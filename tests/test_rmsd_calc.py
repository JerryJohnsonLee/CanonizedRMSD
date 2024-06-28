from CanonizedRMSD import RMSDCalc
import pytest
import json

calc = RMSDCalc()

@pytest.mark.parametrize("name, expected_rmsd, slow", [
    ["a", 6.848138751369438e-16, False],
    ["b", 6.828790187711896e-16, False],
    ["c", 0.08571189099025237, False],
    ["d", 1.6141935076161071, True],
    ["e", 1.6977371875514329, True]
])
def test_rmsd_aligned(name, expected_rmsd, slow, request):
    if not request.config.getoption("--run-slow") and slow:
        pytest.skip(f"Skipping slow test: {name}")
    result = calc.run(f"{pytest.TEST_FOLDER}/{name}/1.sdf", 
                      f"{pytest.TEST_FOLDER}/{name}/2.sdf")
    assert result.rmsd == pytest.approx(expected_rmsd, abs=1e-6)
    with open(f"{pytest.TEST_FOLDER}/{name}/mapping.json") as j:
        ground_truth_mapping = json.load(j)
    assert sorted(result.get_mapping()) == sorted(ground_truth_mapping)

@pytest.mark.parametrize("name, expected_rmsd, slow", [
    ["a", 0, False],
    ["b", 1, False],
    ["c", 0.0856983080346397, False],
    ["d", 4.803715438952769, True],
    ["e", 2.570109531729095, True]
])
def test_rmsd_not_aligned(name, expected_rmsd, slow, request):
    if not request.config.getoption("--run-slow") and slow:
        pytest.skip(f"Skipping slow test: {name}")
    result = calc.run(f"{pytest.TEST_FOLDER}/{name}/1.sdf", 
                      f"{pytest.TEST_FOLDER}/{name}/2.sdf",
                      no_alignment=True)
    assert result.rmsd == pytest.approx(expected_rmsd, abs=1e-6)

