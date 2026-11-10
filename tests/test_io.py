from pathlib import Path
import tempfile
import unittest

import numpy as np

TEMP_ROOT = Path(__file__).resolve().parents[1] / 'tmp/tests/python'
TEMP_ROOT.mkdir(parents=True, exist_ok=True)

from figure.figtools.io import derive_title, discover_energy_series, load_result, parse_multiblock


class ParseMultiblockTests(unittest.TestCase):
    def test_parses_numeric_scientific_notation_and_multiline_headers(self):
        text = """ TOTAL CROSS SECTION AS A FUNCTION OF Z:
 1 8.5529843E-03
 2 4.5

 TOTAL CROSS SECTION AS A FUNCTION OF N:
 1 2.0

 EXCITATION ENERGY CROSS SECTION

 MeV mb/MeV
 0 6.0E-01

 CROSS SECTION OF SPECIAL ELEMETNS:
 8.5529843E-03
"""
        result = parse_multiblock(text)

        np.testing.assert_array_equal(result["Z"][0], [1.0, 2.0])
        np.testing.assert_allclose(result["Z"][1], [8.5529843e-3, 4.5])
        np.testing.assert_allclose(result["excitation_energy"][1], [0.6])
        np.testing.assert_allclose(result["special"][0], [0.0])
        np.testing.assert_allclose(result["special"][1], [8.5529843e-3])

    def test_loads_synthetic_cold_file_and_derives_title(self):
        path = Path("tests/fixtures/CNPCS_COLD_synthetic.dat")
        result = load_result(path)

        self.assertEqual(set(("Z", "N", "A", "special")), set(result))
        self.assertEqual(result["Z"][0][0], 1.0)
        self.assertAlmostEqual(result["Z"][1][1], 4.5)
        self.assertEqual(derive_title(path), "COLD")

    def test_accepts_string_path(self):
        result = parse_multiblock("tests/fixtures/CNPCS_COLD_synthetic.dat")
        self.assertEqual(result.source, Path("tests/fixtures/CNPCS_COLD_synthetic.dat"))

    def test_discovers_energy_series_in_order_and_filters_mode(self):
        with tempfile.TemporaryDirectory(dir=TEMP_ROOT, prefix='io-') as temp_dir:
            root = Path(temp_dir) / "p+232Th"
            for energy in (1000, 80, 300):
                result_dir = root / f"{energy}MeV" / "ana_results"
                result_dir.mkdir(parents=True)
                (result_dir / "CNPCS_COLD.dat").write_text(
                    "TOTAL CROSS SECTION AS A FUNCTION OF Z:\n1 1.0\n",
                    encoding="utf-8",
                )
            hot_dir = root / "200 MeV" / "ana_results"
            hot_dir.mkdir(parents=True)
            (hot_dir / "CNPCS_HOT.dat").write_text(
                "TOTAL CROSS SECTION AS A FUNCTION OF Z:\n1 2.0\n",
                encoding="utf-8",
            )

            cold = discover_energy_series(root, mode="cold")
            hot = discover_energy_series(root, mode="HOT")
            self.assertEqual([label for label, _ in cold], ["80 MeV", "300 MeV", "1000 MeV"])
            self.assertEqual([label for label, _ in hot], ["200 MeV"])
            self.assertEqual(discover_energy_series(root / "missing"), [])


if __name__ == "__main__":
    unittest.main()
