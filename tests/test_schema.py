"""Schema generation is checked against actual parser and serializer behavior."""

import json
from pathlib import Path

import jsonschema
import pytest
from cave_closure import prepare_json

ROOT = Path(__file__).parent.parent


def test_input_result_schemas_and_float_roundtrip():
    problem = json.loads((ROOT / "tests/fixtures/controlled_tree.json").read_text())
    problem["observations"][0]["value"][1] = 0.10000000000000002
    jsonschema.validate(problem, json.loads((ROOT / "schemas/input.json").read_text()))
    result = json.loads(prepare_json(json.dumps(problem)).solve().to_json())
    jsonschema.validate(result, json.loads((ROOT / "schemas/result.json").read_text()))
    assert (
        prepare_json(json.dumps(problem)).solve().to_json()
        == prepare_json(json.dumps(problem)).solve().to_json()
    )


@pytest.mark.parametrize(
    "mutation",
    [
        lambda p: p.update(schema_version=2),
        lambda p: p["frame"].update(something_unknown=True),
        lambda p: p["stations"].__setitem__(0, ["A"]),
        lambda p: p["observations"][0]["uncertainty"].update(undocumented_floor=0.01),
    ],
)
def test_schema_and_parser_both_reject_structural_invalidity(mutation):
    problem = json.loads((ROOT / "tests/fixtures/controlled_tree.json").read_text())
    mutation(problem)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(
            problem, json.loads((ROOT / "schemas/input.json").read_text())
        )
    with pytest.raises(ValueError):
        prepare_json(json.dumps(problem))
