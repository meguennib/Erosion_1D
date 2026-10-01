def test_slfv_policy_is_not_cfl_limited():
    text = open("src/simulation.py", encoding="utf-8").read()
    assert "dt_src, dt_morph, p.dt_max, t_end - t" in text


def test_morphological_accuracy_parameter_is_configured():
    text = open("src/config.py", encoding="utf-8").read()
    assert "morph_rel_change: float = 0.01" in text
