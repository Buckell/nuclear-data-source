#include <iomanip>
#include <iostream>

#include <Python.h>

#include "nds.hpp"

constexpr std::string_view help_message = R"(
=======================================

 ░███    ░██   ░███████      ░██████
 ░████   ░██   ░██   ░██    ░██   ░██
 ░██░██  ░██   ░██    ░██  ░██
 ░██ ░██ ░██   ░██    ░██   ░████████
 ░██  ░██░██   ░██    ░██          ░██
 ░██   ░████   ░██   ░██    ░██   ░██
 ░██    ░███   ░███████      ░██████

         Nuclear Data Source
               v0.0.2
=======================================

  nds - Displays this message.

  Nuclide Format: <Symbol>[Isotope] or <Atomic Number><Isotope> (MCNP format)
  Examples: Co60, Fe, 8016

  nds n <nuclide> - List information for an isotope.

  nds e <expression> - Evaluate a Python expression in the NDS context.
  nds eval <expression>

  nds exec <file> - Execute a Python file in the NDS context.
)";

nds::data_manager dm;

PyObject * Py_nds_nuclide_data(PyObject * self, PyObject * args) {
    const char * nuclide_string;

    if (!PyArg_ParseTuple(args, "s", &nuclide_string)) {
        return nullptr;
    }

    auto const [e, i] = dm.parse_nuclide(nuclide_string);

    auto const & periodic_entry = dm.get_periodic_entry(e);

    PyObject * data = PyDict_New();

    PyDict_SetItemString(data, "name", PyUnicode_FromStringAndSize(periodic_entry.name.data(), static_cast<Py_ssize_t>(periodic_entry.name.size())));

    return data;
}

PyMethodDef Py_nds_methods[] = {
    {"N", Py_nds_nuclide_data, METH_VARARGS, "Get nuclide data." },
    {nullptr, nullptr, 0, nullptr}
};

PyModuleDef Py_nds_module = {
    PyModuleDef_HEAD_INIT, "nds", nullptr, -1, Py_nds_methods
};

PyMODINIT_FUNC PyInit_nds() {
    return PyModule_Create(&Py_nds_module);
}

int main(int const argc, char const * const * const argv) {
    if (argc <= 1) {
        std::cout << help_message << std::endl;
        return 0;
    }

    // auto const endf = ZipFile::Open("./data/ENDF-B-VIII.0_decay.zip");
    // auto const co60 = endf->GetEntry("ENDF-B-VIII.0_decay/dec-027_Co_060.endf");
    // auto const decompress_stream = co60->GetDecompressionStream();
    //
    // std::string line;
    // std::getline(*decompress_stream, line);
    //
    // std::cout << line << '\n';
    //
    // co60->CloseDecompressionStream();

    // Increase double output precision.
    std::cout << std::setprecision(10);

    std::vector<std::string_view> args(argc);
    std::ranges::generate(args, [i = 0, argv]() mutable { return std::string_view(argv[i++]); });

    if (args[1] == "n") {
        auto [e, i] = dm.parse_nuclide(args[2]);

        auto const & periodic_data = dm.get_periodic_entry(e);

        if (i > 0) {
            auto const & nubase_data = dm.get_nubase_entry(e, i);
            auto const awr = dm.get_atomic_weight_ratio(e, i);

            std::cout << "========== " <<  periodic_data.name << '-' << i << " ==========" << '\n';
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
            auto const & nubase_data = dm.get_nubase_entries(e);
            std::ranges::for_each(nubase_data, [](std::pair<std::uint16_t, nds::nubase_entry> const & kvp) {
                auto const & [isotope, data] = kvp;
                std::cout << isotope;

                if (!data.isotope_abundance.empty()) {
                    std::cout << " (" << data.isotope_abundance << "%)";
                }

                std::cout << ", ";
            });
            std::cout << '\n';

            std::cout << "Description: \n" << periodic_data.summary << '\n';

            std::cout << std::flush;

        }
    }
    else if (args[1] == "decay") {
        auto [e, i] = dm.parse_nuclide(args[2]);

        auto const & periodic_data = dm.get_periodic_entry(e);
        auto const decay_data = dm.fetch_decay_data(e, i);

        std::cout << "========== " <<  periodic_data.name << '-' << i << " DECAY DATA ==========" << '\n';

        for (const auto &[energy, energy_delta, intensity, intensity_delta] : decay_data.gamma_discrete) {
            std::cout << nds::format_metric(energy, "eV") << " (" << nds::format_metric(energy_delta, "eV") << ") @ " << intensity * 100  << "% (" << intensity_delta * 100 << "%)\n";
        }

        std::cout << std::flush;
    }
    else if (args[1] == "eval") {
        PyImport_AppendInittab("nds", PyInit_nds);

        Py_Initialize();

        PyObject * main_module = PyImport_AddModule("__main__");
        PyObject * globals = PyModule_GetDict(main_module);
        PyObject * locals = PyDict_New();

        PyObject * func_nuclide_data = PyCFunction_New(Py_nds_methods, nullptr);
        PyDict_SetItemString(globals, "N", func_nuclide_data);

        PyObject * result = PyRun_String("N('Co')['name']", Py_eval_input, globals, locals);

        PyObject * str = PyObject_Str(result);

        char const * buffer = PyUnicode_AsUTF8(str);

        std::cout << buffer << std::endl;

        Py_DECREF(str);
        if (result) Py_DECREF(result);
        Py_Finalize();
    }
}
