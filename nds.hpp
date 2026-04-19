//
// Created by maxng on 12/04/2026.
//

#ifndef NUCLEAR_DATA_SOURCE_NDS_HPP
#define NUCLEAR_DATA_SOURCE_NDS_HPP

#include <algorithm>
#include <array>
#include <charconv>
#include <format>
#include <fstream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ZipFile.h>

constexpr std::string_view periodic_dir = "data/periodic.csv";
constexpr std::size_t periodic_header_offset = 338;

constexpr std::string_view nubase_dir = "data/nubase2020";
constexpr std::size_t nubase_header_offset = 2406;

constexpr std::string_view awr_dir = "data/awr_data";
constexpr std::size_t awr_header_offset = 268;

namespace nds {
    using namespace std::literals;

    constexpr std::double_t neutron_mass_kg = 1.67492750056E-27;
    constexpr std::double_t neutron_mass_mevc2 = 939.56542194;
    constexpr std::double_t neutron_mass_da = 1.008664916;

    inline std::map<std::string_view, std::double_t> nubase_decay_time_unit_factors{
        { "Yy"sv, 1E18 * 365 * 24 * 3600 },
        { "Ey"sv, 1E15 * 365 * 24 * 3600 },
        { "Py"sv, 1E12 * 365 * 24 * 3600 },
        { "Gy"sv, 1E9 * 365 * 24 * 3600 },
        { "My"sv, 1E6 * 365 * 24 * 3600 },
        { "ky"sv, 1E3 * 365 * 24 * 3600 },
        { "y"sv, 365 * 24 * 3600 },
        { "d"sv, 24 * 3600 },
        { "h"sv, 3600 },
        { "m"sv, 60 },
        { "s"sv, 1 },
        { "ms"sv, 1E-3 },
        { "us"sv, 1E-6 },
        { "ns"sv, 1E-9 },
        { "ps"sv, 1E-12 },
        { "fs"sv, 1E-15 },
        { "as"sv, 1E-18 },
        { "zs"sv, 1E-21 },
        { "ys"sv, 1E-24 },
    };

    inline std::map<std::int64_t, std::string_view> metric_prefixes = {
        { 6,  "Y"sv },
        { 5,  "E"sv },
        { 4,  "P"sv },
        { 3,   "G"sv },
        { 2,   "M"sv },
        { 1,   "k"sv },
        { -1,  "m"sv },
        { -2,  "u"sv },
        { -3,  "n"sv },
        { -4, "p"sv },
        { -5, "f"sv },
        { -6, "a"sv },
        { -7, "z"sv },
        { -8, "y"sv },
    };

    constexpr std::string format_metric(std::double_t a_number, std::string_view const a_base_unit) {
        std::int64_t factor = 0;

        while (true) {
            if (a_number > 1000) {
                a_number /= 1000;
                factor += 1;
            } else if (a_number < 0) {
                a_number *= 1000;
                factor -= 1;
            } else {
                break;
            }
        }

        return std::to_string(a_number) + " " + std::string(metric_prefixes[factor]) + std::string(a_base_unit);
    }

    constexpr std::string_view trim(std::string_view a_string) {
        // ReSharper disable once CppDFALocalValueEscapesFunction
        while (!a_string.empty() && a_string[0] == ' ') a_string = a_string.substr(1);
        while (!a_string.empty() && *(a_string.cend() - 1) == ' ') a_string = std::string_view(a_string.cbegin(), a_string.cend() - 1);
        // ReSharper disable once CppDFALocalValueEscapesFunction
        return a_string;
    }

    template <std::integral t_integer>
    constexpr std::optional<t_integer> to_integer(std::string_view const a_string) {
        if (t_integer result; std::from_chars(a_string.data(), a_string.data() + a_string.size(), result).ec == std::errc()) {
            return result;
        }

        return std::nullopt;
    }

    template <std::floating_point t_floating>
    constexpr std::optional<t_floating> to_floating(std::string_view const a_string) {
        if (t_floating result; std::from_chars(a_string.data(), a_string.data() + a_string.size(), result).ec == std::errc()) {
            return result;
        }

        return std::nullopt;
    }

    struct nuclide {
        std::uint8_t const atomic_number;
        std::uint16_t const isotope = 0;
        std::uint8_t const state = 0; // Metastable states, 1 = m, 2 = n
    };

    struct periodic_entry {
        union {
            struct {
                std::string_view const name;
                std::string_view const appearance;
                std::string_view const atomic_mass; // Dalton
                std::string_view const boiling_point; // Kelvin
                std::string_view const category;
                std::string_view const color;
                std::string_view const density; // gas: g/L, solid: g/cc
                std::string_view const discovered_by;
                std::string_view const melting_point; // Kelvin
                std::string_view const molar_heat; // Thermal capacity (J/molK)
                std::string_view const named_by;
                std::string_view const atomic_number;
                std::string_view const period;
                std::string_view const phase;
                std::string_view const source;
                std::string_view const spectral_img; // Link to image.
                std::string_view const summary;
                std::string_view const symbol;
                std::string_view const table_x;
                std::string_view const table_y;
                std::string_view const shells;
                std::string_view const electron_configuration;
                std::string_view const electron_configuration_semantic;
                std::string_view const electron_affinity;
                std::string_view const electronegativity_pauling;
            };

            std::array<std::string_view, 25> fields;
        };

        std::vector<std::string_view> const ionization_energies;

        periodic_entry(std::array<std::string_view, 25> const & a_fields, std::vector<std::string_view> a_ionization_energies) :
            fields(a_fields),
            ionization_energies(std::move(a_ionization_energies))
        {}
    };

    struct nubase_entry {
        std::string_view mass_excess; // keV
        std::string_view half_life;
        std::string_view half_life_units;
        std::string_view isotope_abundance;
        std::double_t decay_constant; // s-1
    };

    struct discrete_energy {
        std::double_t const energy;
        std::double_t const energy_delta;
        std::double_t const intensity;
        std::double_t const intensity_delta;
    };

    struct decay_data {
        std::vector<discrete_energy> const gamma_discrete;
        std::vector<discrete_energy> const xray_discrete;
        std::vector<discrete_energy> const electron_discrete;
    };

    class data_manager {
        std::string periodic_text;
        std::string nubase_text;
        std::string awr_text; // atomic weight ratios

        std::vector<periodic_entry> periodic_data;
        std::vector<std::map<std::uint16_t, std::vector<nubase_entry>>> nubase_data;
        std::vector<std::map<std::uint16_t, std::double_t>> awr_data;

    public:
        data_manager() :
            periodic_text(read_file(periodic_dir)),
            nubase_text(read_file(nubase_dir)),
            awr_text(read_file(awr_dir)),
            nubase_data(119),
            awr_data(119)
        {
            periodic_data.reserve(119);

            // Database cleanup and processing.
            process_periodic();
            process_nubase();
            process_awr();
        }

        [[nodiscard]]
        std::double_t get_atomic_weight_ratio(nuclide const a_nuclide) const {
            return awr_data[a_nuclide.atomic_number - 1].at(a_nuclide.isotope);
        }

        [[nodiscard]]
        periodic_entry const & get_periodic_entry(std::uint8_t const a_atomic_number) const {
            return periodic_data[a_atomic_number - 1];
        }

        [[nodiscard]]
        std::map<std::uint16_t, std::vector<nubase_entry>> const & get_nubase_entries(std::uint8_t const a_atomic_number) const {
            return nubase_data[a_atomic_number - 1];
        }

        [[nodiscard]]
        nubase_entry const & get_nubase_entry(nuclide const a_nuclide) const {
            auto const [atomic_number, isotope, state] = a_nuclide;
            return nubase_data[atomic_number - 1].at(isotope).at(state);
        }

        [[nodiscard]]
        nuclide parse_nuclide(std::string_view a_nuclide) const {
            auto state = 0ui8;

            if (a_nuclide.ends_with('m')) {
                state = 1;
                a_nuclide = a_nuclide.substr(0, a_nuclide.length() - 1);
            }
            else if (a_nuclide.ends_with('n')) {
                state = 2;
                a_nuclide = a_nuclide.substr(0, a_nuclide.length() - 1);
            }

            auto isotope_start = 0;
            while (isotope_start < a_nuclide.length() && !(a_nuclide[isotope_start] >= '0' && a_nuclide[isotope_start] <= '9')) ++isotope_start;

            auto element_string = "000"sv;
            auto isotope_string = "000"sv;

            if (isotope_start == 0) {
                isotope_string = a_nuclide.substr(a_nuclide.length() - 3);
                element_string = a_nuclide.substr(0, a_nuclide.length() - isotope_string.length());
            } else {
                auto element_symbol = a_nuclide;
                isotope_string = "000"sv;

                if (isotope_start != a_nuclide.length()) {
                    element_symbol = a_nuclide.substr(0, isotope_start);
                    isotope_string = a_nuclide.substr(isotope_start);
                }

                auto const periodic_it = std::ranges::find_if(periodic_data, [element_symbol](periodic_entry const & a_entry) {
                    return a_entry.symbol == element_symbol;
                });

                if (periodic_it == periodic_data.cend()) {
                    return {0, 0};
                }

                element_string = periodic_it->atomic_number;
            }

            std::uint8_t element;
            std::uint16_t isotope;

            std::from_chars(element_string.data(), element_string.data() + element_string.size(), element);
            std::from_chars(isotope_string.data(), isotope_string.data() + isotope_string.size(), isotope);

            return { element, isotope, state };
        }

        [[nodiscard]]
        decay_data fetch_decay_data(nuclide const a_nuclide) const {
            auto const [atomic_number, isotope, state] = a_nuclide;

            auto const convert_endf_numeral = [](std::string_view const a_numeral_string) -> std::double_t {
                if (auto const plus_position = a_numeral_string.rfind('+'); plus_position != std::string_view::npos) {
                    return to_floating<std::double_t>(a_numeral_string.substr(0, plus_position)).value()
                           * pow(10, to_integer<std::uint8_t>(a_numeral_string.substr(plus_position + 1)).value());
                }

                if (auto const negative_position = a_numeral_string.rfind('-'); negative_position != std::string_view::npos) {
                    return to_floating<std::double_t>(a_numeral_string.substr(0, negative_position)).value()
                           * pow(10, -to_integer<std::uint8_t>(a_numeral_string.substr(negative_position + 1)).value());
                }

                return 0;
            };

            auto const & periodic_entry = get_periodic_entry(atomic_number);

            auto atomic_number_string = std::to_string(atomic_number);
            atomic_number_string = std::string(3 - atomic_number_string.length(), '0') + atomic_number_string;
            auto isotope_string = std::to_string(isotope);
            isotope_string = std::string(3 - isotope_string.length(), '0') + isotope_string;
            auto state_string = state > 0 ? std::format("m{}", state) : "";

            auto const endf = ZipFile::Open("./data/ENDF-B-VIII.0_decay.zip");
            auto const decay_archive = endf->GetEntry(std::format("ENDF-B-VIII.0_decay/dec-{}_{}_{}{}.endf", atomic_number_string, periodic_entry.symbol, isotope_string, state_string));
            auto const decompress_stream = decay_archive->GetDecompressionStream();

            std::vector<discrete_energy> gammas;
            std::vector<discrete_energy> xrays;
            std::vector<discrete_energy> electrons;

            for (std::string line; !(line.length() > 68 && line[68] == '-'); std::getline(*decompress_stream, line)) {
                if (std::string_view line_view = line; line_view.length() > 71 && line_view.substr(71, 4) == "8457") {
                    if (
                        auto const terminal_value_string = trim(line_view.substr(56, 10));
                        terminal_value_string != "0" && terminal_value_string.rfind('-') == std::string_view::npos && terminal_value_string.rfind('+') == std::string_view::npos
                    ) {
                        // 0: gamma
                        // 1: B- endpoints
                        // 8: electron
                        // 9: X-ray
                        auto const key = convert_endf_numeral(line.substr(12, 10));

                        if (convert_endf_numeral(line.substr(1, 10)) != 0) {
                            continue;
                        }

                        auto section_entry_count = to_integer<std::size_t>(trim(line_view.substr(56, 10))).value();

                        for (auto rows_left = to_integer<std::size_t>(trim(line_view.substr(45, 10))).value() / 6; rows_left > 0; --rows_left) {
                            std::getline(*decompress_stream, line);
                        }

                        while (section_entry_count--) {
                            std::getline(*decompress_stream, line);
                            line_view = line;

                            auto rows_left = to_integer<std::size_t>(trim(line_view.substr(45, 10))).value() / 6;

                            auto const energy = convert_endf_numeral(line_view.substr(1, 10));
                            auto const energy_delta = convert_endf_numeral(line_view.substr(12, 10));

                            std::getline(*decompress_stream, line);
                            line_view = line;

                            --rows_left;

                            auto const intensity = convert_endf_numeral(line_view.substr(23, 10));
                            auto const intensity_delta = convert_endf_numeral(line_view.substr(34, 10));

                            switch (static_cast<std::size_t>(key)) {
                                case 0:
                                    gammas.emplace_back(energy, energy_delta, intensity, intensity_delta);
                                    break;
                                case 8:
                                    electrons.emplace_back(energy, energy_delta, intensity, intensity_delta);
                                    break;
                                case 9:
                                    xrays.emplace_back(energy, energy_delta, intensity, intensity_delta);
                                    break;
                                default:
                                    break;
                            }

                            for (; rows_left > 0; --rows_left) {
                                std::getline(*decompress_stream, line);
                            }
                        }
                    }
                }
            }

            decay_archive->CloseDecompressionStream();

            return {
                .gamma_discrete = gammas,
                .xray_discrete = xrays,
                .electron_discrete = electrons,
            };
        }

    private:
        void process_periodic() {
            // Trim header.
            periodic_text.erase(periodic_text.cbegin(), periodic_text.cbegin() + periodic_header_offset);

            // Quickly process CSV.
            std::iter_difference_t<std::string_view> scan_index = 0;
            std::array<std::iter_difference_t<std::string_view>, 25> comma_indices{};

            // Scan all element lines.
            for (std::size_t i = 0; i < 119; ++i) {
                std::size_t current_comma_index = 0;

                auto const line_start = scan_index;

                // Locate commas.

                bool quote_mode = false;
                for (auto scan_char = periodic_text[scan_index]; scan_char != '\n'; scan_char = periodic_text[++scan_index]) {
                    if (scan_char == '"') {
                        quote_mode = !quote_mode;
                    } else if (scan_char == ',' && !quote_mode) {
                        comma_indices[current_comma_index++] = scan_index;
                    }

                    if (scan_index == periodic_text.size()) {
                        break;
                    }
                }
                ++scan_index;

                // Split into strings.

                std::array<std::string_view, 25> fields{std::string_view(periodic_text.cbegin() + line_start + 1, periodic_text.cbegin() + comma_indices[0] - 1)};
                std::generate(fields.begin() + 1, fields.end() - 1, [this, &comma_indices, ii = 0]() mutable {
                    auto ret = std::string_view(periodic_text.cbegin() + comma_indices[ii] + 1, periodic_text.cbegin() + comma_indices[ii + 1]); // NOLINT(*-dangling-handle)
                    if (ret.starts_with('"')) {
                        ret = ret.substr(1, ret.size() - 2);
                    }

                    ++ii;
                    // ReSharper disable once CppDFALocalValueEscapesFunction
                    return ret;
                });
                *(fields.end() - 1) = std::string_view(periodic_text.cbegin() + comma_indices[23] + 1, periodic_text.cbegin() + comma_indices[24]); // NOLINT(*-dangling-handle)

                // Find ionization energies.

                std::vector<std::string_view> ionization_energies;

                auto base = periodic_text.cbegin() + *(comma_indices.cend() - 1) + 3;

                if (!std::string_view(base, periodic_text.cbegin() + (scan_index - 3)).empty()) {
                    for (auto it = base; it != periodic_text.cend() && *it != ']'; ++it) {
                        if (*it == ',') {
                            ionization_energies.emplace_back(base, it);
                            base = it + 1;
                        }
                    }
                    ionization_energies.emplace_back(base, periodic_text.cbegin() + (scan_index - 3));
                }

                periodic_data.emplace_back(fields, ionization_energies);
            }
        }

        void process_nubase() {
            // Trim heading.
            nubase_text.erase(nubase_text.cbegin(), nubase_text.cbegin() + nubase_header_offset);

            auto const process_line = [this](std::string_view const a_line) {
                auto const mass_number   = a_line.substr(0, 3);
                auto const atomic_number = a_line.substr(4, 3);
                auto const alt_type      = a_line.substr(8, 1);
                auto const isotope_name  = a_line.substr(11, 5);
                auto const alt           = a_line.substr(16, 1);
                auto const mass_excess   = trim(a_line.substr(18, 13)); // keV

                auto const half_life = a_line.length() > 69 ? trim(a_line.substr(69, 9)) : ""sv;
                auto const half_life_units = a_line.length() > 78 ? trim(a_line.substr(78, 2)) : ""sv;

                // Scan decay/isotope extra information.

                auto const decay_extra_info = a_line.length() > 119 ? a_line.substr(119) : ""sv;

                auto isotope_abundance = ""sv;

                auto decay_ei_entry_start = decay_extra_info.cbegin();
                for (auto it = decay_ei_entry_start; it <= decay_extra_info.cend(); ++it) {
                    if (it == decay_extra_info.cend() || *it == ';') {
                        auto const entry = std::string_view(decay_ei_entry_start, it);

                        if (entry.starts_with("IS=")) {
                            isotope_abundance = entry.substr("IS="sv.length());
                        }

                        if (it == decay_extra_info.cend()) {
                            break;
                        }

                        decay_ei_entry_start = it;
                    }
                }

                auto const element = to_integer<std::uint8_t>(atomic_number).value();
                auto const isotope = to_integer<std::uint16_t>(mass_number).value();
                auto const state = to_integer<std::uint16_t>(alt_type).value_or(0);

                if (element == 0) {
                    return;
                }

                auto const half_life_value = to_floating<std::double_t>(half_life);

                auto & element_nubase_entries = nubase_data[element - 1];
                auto [pair_it, did_emplace] = element_nubase_entries.try_emplace(isotope, std::vector<nubase_entry>());
                auto & [entry_isotope, entry_vector] = *pair_it;

                if (entry_vector.size() <= state) {
                    entry_vector.resize(state);
                }

                entry_vector.emplace(entry_vector.cbegin() + state,
                    mass_excess,
                    half_life,
                    half_life_units,
                    isotope_abundance,
                    half_life_value.has_value() ? log(2) / (half_life_value.value() * nubase_decay_time_unit_factors[half_life_units]) : 0.0
                );
            };

            auto line_begin = nubase_text.cbegin();
            for (auto it = nubase_text.cbegin(); it != nubase_text.cend(); ++it) {
                if (*it == '\n') {
                    process_line({ line_begin, it });
                    line_begin = it + 1;
                }
            }
        }

        void process_awr() {
            // Trim header.
            awr_text.erase(awr_text.cbegin(), awr_text.cbegin() + awr_header_offset);

            std::vector<std::string_view> entry_lines;

            auto scan_it = awr_text.cbegin();
            auto line_begin = scan_it;
            for (; scan_it != awr_text.cend(); ++scan_it) {
                if (*scan_it == '\n') {
                    auto const current_line = std::string_view(line_begin, scan_it); // NOLINT(*-dangling-handle)
                    line_begin = scan_it + 1;

                    if (current_line.substr(0, 7) != "       ") {
                        if (!entry_lines.empty()) {
                            // Process entry.

                            auto const atomic_number_string = trim(entry_lines[0].substr(1, 3));
                            auto const atomic_number = to_integer<std::uint8_t>(atomic_number_string);

                            if (!atomic_number.has_value() || atomic_number == 0) {
                                entry_lines.clear();
                                entry_lines.emplace_back(current_line);
                                continue;
                            }

                            auto & element_data = awr_data[atomic_number.value() - 1];

                            element_data[0] = to_floating<std::double_t>(entry_lines[0].substr(9, 10)).value();

                            for (auto const line : entry_lines) {
                                if (line.length() > 27) {
                                    auto const first_isotope = to_integer<std::uint16_t>(trim(line.substr(20 + atomic_number_string.length(), 3))).value();
                                    if (line.length() > 47) {
                                        auto const second_isotope = to_integer<std::uint16_t>(trim(line.substr(40 + atomic_number_string.length(), 3))).value();

                                        if (line.length() > 67) {
                                            auto const third_isotope = to_integer<std::uint16_t>(trim(line.substr(60 + atomic_number_string.length(), 3))).value();
                                            auto const third_ratio = to_floating<std::double_t>(trim(line.substr(67))).value();
                                            element_data[third_isotope] = third_ratio;

                                            auto const second_ratio = to_floating<std::double_t>(trim(line.substr(47, 12))).value();
                                            element_data[second_isotope] = second_ratio;
                                        } else {
                                            auto const second_ratio = to_floating<std::double_t>(trim(line.substr(47))).value();
                                            element_data[second_isotope] = second_ratio;
                                        }

                                        auto const first_ratio = to_floating<std::double_t>(trim(line.substr(27, 12))).value();
                                        element_data[first_isotope] = first_ratio;
                                    } else {
                                        auto const first_ratio = to_floating<std::double_t>(trim(line.substr(27))).value();
                                        element_data[first_isotope] = first_ratio;
                                    }
                                }
                            }

                            entry_lines.clear();
                        }
                    }

                    entry_lines.emplace_back(current_line);
                }
            }
        }

        [[nodiscard]]
        static std::string read_file(std::string_view const a_directory) {
            return read_file(std::string(a_directory));
        }

        [[nodiscard]]
        static std::string read_file(std::string const & a_directory) {
            std::ifstream const file_stream(a_directory);
            std::stringstream string_stream;
            string_stream << file_stream.rdbuf();

            return string_stream.str();
        }
    };
}

#endif //NUCLEAR_DATA_SOURCE_NDS_HPP
