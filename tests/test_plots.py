from pathlib import Path
import tempfile
import unittest

TEMP_ROOT = Path(__file__).resolve().parents[1] / 'tmp/tests/python'
TEMP_ROOT.mkdir(parents=True, exist_ok=True)

from matplotlib.figure import Figure

from figure.figtools.io import Result
from figure.figtools.plots import plot_cross_section_distributions, plot_energy_comparison


class PlotTests(unittest.TestCase):
    def test_returns_independent_figures_and_hides_missing_blocks(self):
        result = Result({"Z": ([1.0, 2.0], [3.0, 4.0]), "A": ([10.0], [5.0])})
        with tempfile.TemporaryDirectory(dir=TEMP_ROOT, prefix='plots-') as temp_dir:
            first_path = Path(temp_dir) / "first.png"
            second_path = Path(temp_dir) / "second.png"
            first = plot_cross_section_distributions(result, output=first_path)
            second = plot_cross_section_distributions(result, output=second_path, quantities=("A",), yscale="log")

            self.assertIsInstance(first, Figure)
            self.assertIsInstance(second, Figure)
            self.assertIsNot(first, second)
            self.assertEqual(len(first.axes), 3)
            self.assertFalse(first.axes[1].get_visible())
            self.assertTrue(first_path.exists())
            self.assertTrue(second_path.exists())
            self.assertEqual(second.axes[0].get_yscale(), "log")

    def test_energy_comparison_has_six_plain_lines_and_log_axis(self):
        series = [
            (f"{energy} MeV", Result({"A": ([1.0, 2.0], [energy, energy / 2])}))
            for energy in (80, 100, 200, 300, 600, 1000)
        ]
        with tempfile.TemporaryDirectory(dir=TEMP_ROOT, prefix='plots-') as temp_dir:
            output = Path(temp_dir) / "A_energy.png"
            figure = plot_energy_comparison(series, "A", output=output, title="A comparison")

            self.assertIsInstance(figure, Figure)
            self.assertEqual(len(figure.axes), 1)
            self.assertEqual(figure.axes[0].get_yscale(), "log")
            self.assertEqual(len(figure.axes[0].lines), 6)
            self.assertTrue(all(line.get_marker() in ("None", None) for line in figure.axes[0].lines))
            self.assertEqual(len(figure.axes[0].get_legend().get_texts()), 6)
            self.assertTrue(output.exists())


if __name__ == "__main__":
    unittest.main()
