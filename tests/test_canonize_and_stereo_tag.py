from CanonizedRMSD import Canonizer
import pytest
import json

canonizer = Canonizer()
with open(f"{pytest.TEST_FOLDER}/ground_truth.json") as f:
    ground_truth = json.load(f)

@pytest.mark.parametrize("name", ["g", "h", "i", "j", "k", "l"])
def test_rmsd_aligned(name):
    result = canonizer.canonize(f"{pytest.TEST_FOLDER}/{name}.mol")
    assert result.original_to_canonized_mapping == ground_truth[name]["original_to_canonized_mapping"]
    
    stereo_tags = result.get_stereochemical_tags()
    str_named_tags = {str(k): v for k, v in stereo_tags.items()}
    assert str_named_tags == ground_truth[name]["stereo_tags"]