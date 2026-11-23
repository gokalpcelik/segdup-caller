import pytest
from genecaller import genecaller as gc
from genecaller.util import merge_gene_vcfs


class StubGene:
    def __init__(self, name, fail=False):
        self.gene_names = [name]
        self.config = {}
        self.fail = fail


def stub_process(gene_data, gene_logger):
    gene = gene_data[0]
    if gene.fail:
        raise RuntimeError("boom")
    return {"Copy numbers": {gene.gene_names[0]: 2}}


@pytest.fixture(autouse=True)
def patch(monkeypatch):
    monkeypatch.setattr(gc, "_process_single_gene", stub_process)
    monkeypatch.setattr(gc, "_source_cnvscope_stage", lambda *a, **k: {})


def make_caller(genes, workers):
    c = gc.GeneCaller.__new__(gc.GeneCaller)
    c.genes = genes
    c.gene_sections = [g.gene_names[0] for g in genes]
    c.input = None
    c.ref = None
    c.params = {"workers": workers, "threads": 1}
    return c


def test_failed_gene_result_shape():
    r = gc._failed_gene_result("GBA", ValueError("bad"))
    assert r["status"] == "failed"
    assert "GBA" in r["note"] and "ValueError" in r["note"]
    assert r["error"]["type"] == "ValueError"


@pytest.mark.parametrize("workers", [1, 2])
def test_one_failure_does_not_stop_others(workers):
    genes = [StubGene("A"), StubGene("B", fail=True), StubGene("C")]
    res = make_caller(genes, workers).process_genes()
    assert list(res) == ["A", "B", "C"]
    assert res["B"]["status"] == "failed" and "note" in res["B"]
    assert res["A"]["Copy numbers"] == {"A": 2}
    assert res["C"]["Copy numbers"] == {"C": 2}


def test_merge_ignores_failed_entry(tmp_path):
    with pytest.raises(FileNotFoundError):
        merge_gene_vcfs({"B": gc._failed_gene_result("B", RuntimeError("x"))}, str(tmp_path / "m.vcf"))
