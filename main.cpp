#include <iomanip>
#include <iostream>
#include <numeric>

#include <Python.h>
#include <ranges>
#include <span>

#include "nds.hpp"

constexpr std::wstring_view intro_message = LR"(
====================================================================

  NNNNNNNN        NNNNNNNNDDDDDDDDDDDDD           SSSSSSSSSSSSSSS
  N:::::::N       N::::::ND::::::::::::DDD      SS:::::::::::::::S
  N::::::::N      N::::::ND:::::::::::::::DD   S:::::SSSSSS::::::S
  N:::::::::N     N::::::NDDD:::::DDDDD:::::D  S:::::S     SSSSSSS
  N::::::::::N    N::::::N  D:::::D    D:::::D S:::::S
  N:::::::::::N   N::::::N  D:::::D     D:::::DS:::::S
  N:::::::N::::N  N::::::N  D:::::D     D:::::D S::::SSSS
  N::::::N N::::N N::::::N  D:::::D     D:::::D  SS::::::SSSSS
  N::::::N  N::::N:::::::N  D:::::D     D:::::D    SSS::::::::SS
  N::::::N   N:::::::::::N  D:::::D     D:::::D       SSSSSS::::S
  N::::::N    N::::::::::N  D:::::D     D:::::D            S:::::S
  N::::::N     N:::::::::N  D:::::D    D:::::D             S:::::S
  N::::::N      N::::::::NDDD:::::DDDDD:::::D  SSSSSSS     S:::::S
  N::::::N       N:::::::ND:::::::::::::::DD   S::::::SSSSSS:::::S
  N::::::N        N::::::ND::::::::::::DDD     S:::::::::::::::SS
  NNNNNNNN         NNNNNNNDDDDDDDDDDDDD         SSSSSSSSSSSSSSS

                      Nuclear Data Source
                            v0.1.1
====================================================================

    Type 'help' to view available commands.
)";

constexpr std::string_view help_message = R"(
  help - Displays this message.

  n <nuclide> - List information for an isotope.
  decay <nuclide> - List decay information for an isotope (discrete energies).

  Nuclide Format: <Symbol>[Isotope] or <Atomic Number><Isotope> (MCNP format)
  Examples: Co60, Fe, 8016

  eval <expression> - Evaluate a Python expression within the NDS context.

  exit - Leave the program.
)";

//   nds exec <file> - Execute a Python file in the NDS context.

nds::data_manager dm;

bool python_mode = false;

PyObject * main_module;
PyObject * globals;
PyObject * locals;

void execute_command(std::span<std::string_view const> a_args);
std::vector<std::string> parse_arguments(std::string_view a_command);

void initialize_python_environment();
void deinitialize_python_environment();
void run_python_command(std::string const & a_command);

PyObject * Py_nds_nuclide_data(PyObject * self, PyObject * args);
PyObject * Py_nds_nuclide_decay_data(PyObject * self, PyObject * args);

PyMethodDef Py_nds_methods[] = {
    {"N", Py_nds_nuclide_data, METH_VARARGS, "Get nuclide data." },
    {"DecayData", Py_nds_nuclide_decay_data, METH_VARARGS, "Get nuclide decay data." },
    {nullptr, nullptr, 0, nullptr}
};

int main(int const argc, char const * const * const argv) {
    std::vector<std::string_view> initial_args(argc);
    std::ranges::generate(initial_args, [i = 0, argv]() mutable { return std::string_view(argv[i++]); });

    initialize_python_environment();

    // Increase floating point output precision.
    std::cout << std::setprecision(10);

    if (argc > 1) {
        execute_command(std::span(initial_args).subspan(1));
        return EXIT_SUCCESS;
    }

    std::wcout << intro_message << std::endl;

    std::string command;
    do {
        if (python_mode) {
            std::cout << "PY ";
        }

        std::cout << ">>> " << std::flush;
        std::getline(std::cin, command);

        if (command == "py") {
            python_mode = !python_mode;

            if (python_mode) {
                std::cout << "Python mode enabled. Type 'py' again to quit." << std::endl;
            }
            continue;
        }

        if (python_mode) {
            run_python_command(command);
            continue;
        }

        auto const args = parse_arguments(command);
        auto const args_view = args
            | std::views::transform([](auto const & a_arg) { return std::string_view(a_arg); })
            | std::ranges::to<std::vector<std::string_view>>();
        execute_command(args_view);
    } while (command != "exit");

    std::cout << "Thank you for using NDS!" << std::endl;

    deinitialize_python_environment();

    return EXIT_SUCCESS;
}

// ==================================== COMMAND HANDLING ====================================

void execute_command(std::span<std::string_view const> const a_args) {
    if (a_args[0] == "n") {
        auto const nuclide = dm.parse_nuclide(a_args[1]);

        auto const & periodic_data = dm.get_periodic_entry(nuclide.atomic_number);

        if (nuclide.isotope > 0) {
            auto const & nubase_data = dm.get_nubase_entry(nuclide);
            auto const awr = dm.get_atomic_weight_ratio(nuclide);

            std::cout << "========== " <<  periodic_data.name << '-' << nuclide.isotope << " ==========" << '\n';
            std::cout << "Atomic Mass: " << awr * nds::neutron_mass_da << " u (" << awr << "n)\n";
            std::cout << "             " << awr * nds::neutron_mass_mevc2 << " MeV/c^2\n";
            std::cout << "             " << awr * nds::neutron_mass_kg << " kg\n";
            if (!nubase_data.isotope_abundance.empty()) std::cout << "Natural Abundance: " << nubase_data.isotope_abundance << '%' << '\n';
            std::cout << "Mass Excess: " << nubase_data.mass_excess << " keV" << '\n';
            if (!nubase_data.half_life.empty() && nubase_data.half_life != "stbl") {
                std::cout << "Half-Life: " << nubase_data.half_life << ' ' << nubase_data.half_life_units << '\n';
                std::cout << "Decay Constant: " << nubase_data.decay_constant << " s^-1 \n";
            }

            std::cout << std::flush;
        } else {
            std::cout << "========== " <<  periodic_data.name << " (" << periodic_data.atomic_number << ") ==========" << '\n';

            std::cout << "Atomic Mass: " << periodic_data.atomic_mass << " u \n";

            if (periodic_data.phase == "Gas") {
                std::cout << "Density: " << periodic_data.density << " g/L \n";
            } else {
                std::cout << "Density: " << periodic_data.density << " g/cc \n";
            }

            std::cout << "Boiling Point: " << periodic_data.boiling_point << " K \n";
            std::cout << "Melting Point: " << periodic_data.melting_point << " K \n";

            std::cout << "Isotopes: ";
            auto const & nubase_data = dm.get_nubase_entries(nuclide.atomic_number);
            std::ranges::for_each(nubase_data, [](std::pair<std::uint16_t, std::vector<nds::nubase_entry>> const & kvp) {
                auto const & [isotope, data] = kvp;
                std::cout << isotope;

                if (!data[0].isotope_abundance.empty()) {
                    std::cout << " (" << data[0].isotope_abundance << "%)";
                }

                std::cout << ", ";
            });
            std::cout << '\n';

            std::cout << "Description: \n" << periodic_data.summary << '\n';

            std::cout << std::flush;

        }
    }
    else if (a_args[0] == "decay") {
        auto const nuclide = dm.parse_nuclide(a_args[1]);

        auto const & periodic_data = dm.get_periodic_entry(nuclide.atomic_number);
        const auto [gamma_discrete, xray_discrete, electron_discrete] = dm.fetch_decay_data(nuclide);

        std::cout << "========== " <<  periodic_data.name << '-' << nuclide.isotope << " DECAY DATA ==========" << '\n';

        std::cout << "GAMMA: \n";
        for (const auto &[energy, energy_delta, intensity, intensity_delta] : gamma_discrete) {
            std::cout << nds::format_metric(energy, "eV") << " (" << nds::format_metric(energy_delta, "eV") << ") @ " << intensity * 100  << "% (" << intensity_delta * 100 << "%)\n";
        }

        std::cout << "X-RAY: \n";
        for (const auto &[energy, energy_delta, intensity, intensity_delta] : xray_discrete) {
            std::cout << nds::format_metric(energy, "eV") << " (" << nds::format_metric(energy_delta, "eV") << ") @ " << intensity * 100  << "% (" << intensity_delta * 100 << "%)\n";
        }

        std::cout << "ELECTRON: \n";
        for (const auto &[energy, energy_delta, intensity, intensity_delta] : electron_discrete) {
            std::cout << nds::format_metric(energy, "eV") << " (" << nds::format_metric(energy_delta, "eV") << ") @ " << intensity * 100  << "% (" << intensity_delta * 100 << "%)\n";
        }

        std::cout << std::flush;
    }
    else if (a_args[0] == "eval") {
        std::string command;

        std::ranges::for_each(a_args, [&command](auto const & a_arg) { (command += a_arg) += ' '; });

        run_python_command(command.substr(5));
    }
}

std::vector<std::string> parse_arguments(std::string_view const a_command) {
    std::vector<std::string> args;

    bool escape_mode = false;
    bool quote_mode = false;
    auto arg_start_it = a_command.cbegin();

    std::vector<std::size_t> exclusions;

    for (auto it = arg_start_it; it != a_command.cend(); ++it) {
        if (escape_mode) {
            escape_mode = false;
            continue;
        };

        if (*it == '\\') {
            escape_mode = true;
            exclusions.push_back(it - arg_start_it);
            continue;
        }

        if (*it == '"') {
            if (quote_mode) {
                if (it + 1 != a_command.cend() && *(it + 1) == ' ') {
                    quote_mode = false;

                    auto arg = std::string(arg_start_it, it);
                    std::ranges::for_each(exclusions, [i = 0, &arg](auto const exclusion) mutable {
                        arg.erase(exclusion - i++, 1);
                    });
                    exclusions.clear();

                    args.emplace_back(std::move(arg));
                    arg_start_it = ++it + 1;
                }
            } else if (it == arg_start_it) {
                arg_start_it = it + 1;
                quote_mode = true;
            }

            continue;
        }

        if (quote_mode) continue;

        if (*it == ' ') {
            auto arg = std::string(arg_start_it, it);
            std::ranges::for_each(exclusions, [i = 0, &arg](auto const exclusion) mutable {
                arg.erase(exclusion - i++, 1);
            });
            exclusions.clear();

            args.emplace_back(std::move(arg));
            arg_start_it = it + 1;
        }
    }

    if (arg_start_it != a_command.cend()) {
        auto arg = std::string(arg_start_it, a_command.cend());
        std::ranges::for_each(exclusions, [i = 0, &arg](auto const exclusion) mutable {
            arg.erase(exclusion - i++, 1);
        });
        exclusions.clear();

        args.emplace_back(std::move(arg));
    }

    return std::move(args);
}


// ==================================== PYTHON ENVIRONMENT ====================================

void initialize_python_environment() {
    Py_Initialize();

    main_module = PyImport_AddModule("__main__");
    globals = PyModule_GetDict(main_module);
    locals = PyDict_New();

    for (auto ptr = Py_nds_methods; ptr->ml_meth != nullptr; ++ptr) {
        PyObject * func = PyCFunction_New(ptr, nullptr);
        PyDict_SetItemString(globals, ptr->ml_name, func);
        Py_DECREF(func);
    }
}

void deinitialize_python_environment() {
    Py_Finalize();
}

void run_python_command(std::string const & a_command) {
    if (PyObject const * const result = PyRun_String(a_command.data(), Py_single_input, globals, locals); result == nullptr) {
        if (PyErr_Occurred()) {
            PyErr_Print();

            PyObject * type, * value, * traceback;
            PyErr_Fetch(&type, &value, &traceback);

            PyErr_NormalizeException(&type, &value, &traceback);

            Py_XDECREF(type);
            Py_XDECREF(value);
            Py_XDECREF(traceback);
        }

        return;
    }

    // PyObject * str = PyObject_Str(result);
    //
    // char const * buffer = PyUnicode_AsUTF8(str);
    //
    // std::cout << buffer << std::endl;
    //
    // Py_DECREF(str);
    // Py_XDECREF(result);
}

// ==================================== PYTHON METHODS ====================================

PyObject * Py_nds_nuclide_data(PyObject *, PyObject * args) {
    const char * nuclide_string;

    if (!PyArg_ParseTuple(args, "s", &nuclide_string)) {
        return nullptr;
    }

    auto const nuclide = dm.parse_nuclide(nuclide_string);

    auto const & periodic_entry = dm.get_periodic_entry(nuclide.atomic_number);

    PyObject * data = PyDict_New();

    PyDict_SetItemString(data, "name", PyUnicode_FromStringAndSize(periodic_entry.name.data(), static_cast<Py_ssize_t>(periodic_entry.name.size())));

    return data;
}

PyObject * Py_nds_nuclide_decay_data(PyObject * self, PyObject * args) {
    const char * nuclide_string;

    if (!PyArg_ParseTuple(args, "s", &nuclide_string)) {
        return nullptr;
    }

    auto const nuclide = dm.parse_nuclide(nuclide_string);

    const auto [gamma_discrete, xray_discrete, electron_discrete] = dm.fetch_decay_data(nuclide);

    PyObject * py_data = PyDict_New();

    PyObject * py_gammas = PyList_New(static_cast<Py_ssize_t>(gamma_discrete.size()));
    Py_ssize_t py_gamma_ctr = 0;

    for (const auto [energy, energy_delta, intensity, intensity_delta] : gamma_discrete) {
        PyObject * py_row = PyList_New(4);

        PyList_SetItem(py_row, 0, PyFloat_FromDouble(energy));
        PyList_SetItem(py_row, 1, PyFloat_FromDouble(energy_delta));
        PyList_SetItem(py_row, 2, PyFloat_FromDouble(intensity));
        PyList_SetItem(py_row, 3, PyFloat_FromDouble(intensity_delta));

        PyList_SetItem(py_gammas, py_gamma_ctr++, py_row);
    }

    PyObject * py_xrays = PyList_New(static_cast<Py_ssize_t>(xray_discrete.size()));
    Py_ssize_t py_xray_ctr = 0;

    for (const auto [energy, energy_delta, intensity, intensity_delta] : xray_discrete) {
        PyObject * py_row = PyList_New(4);

        PyList_SetItem(py_row, 0, PyFloat_FromDouble(energy));
        PyList_SetItem(py_row, 1, PyFloat_FromDouble(energy_delta));
        PyList_SetItem(py_row, 2, PyFloat_FromDouble(intensity));
        PyList_SetItem(py_row, 3, PyFloat_FromDouble(intensity_delta));

        PyList_SetItem(py_xrays, py_xray_ctr++, py_row);
    }

    PyObject * py_electrons = PyList_New(static_cast<Py_ssize_t>(electron_discrete.size()));
    Py_ssize_t py_electron_ctr = 0;

    for (const auto [energy, energy_delta, intensity, intensity_delta] : electron_discrete) {
        PyObject * py_row = PyList_New(4);

        PyList_SetItem(py_row, 0, PyFloat_FromDouble(energy));
        PyList_SetItem(py_row, 1, PyFloat_FromDouble(energy_delta));
        PyList_SetItem(py_row, 2, PyFloat_FromDouble(intensity));
        PyList_SetItem(py_row, 3, PyFloat_FromDouble(intensity_delta));

        PyList_SetItem(py_electrons, py_electron_ctr++, py_row);
    }

    PyDict_SetItemString(py_data, "gamma_discrete", py_gammas);
    Py_DECREF(py_gammas);

    PyDict_SetItemString(py_data, "xray_discrete", py_xrays);
    Py_DECREF(py_xrays);

    PyDict_SetItemString(py_data, "electron_discrete", py_electrons);
    Py_DECREF(py_electrons);

    return py_data;
}