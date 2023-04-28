import pytest
from paraphase.paraphase import Paraphase


class TestParaphase(object):
    def test_get_gene_list(self):
        app = Paraphase()
        app.accepted_genes = ["smn1", "ncf1"]

        app.genes_to_call = None
        gene_input = "smn1"
        gene_list = app.get_gene_list(gene_input)
        assert gene_list == ["smn1"]

        app.genes_to_call = None
        gene_input = "ncf1,smn1"
        gene_list = app.get_gene_list(gene_input)
        assert sorted(gene_list) == ["ncf1", "smn1"]

        app.genes_to_call = ["smn1"]
        gene_input = None
        gene_list = app.get_gene_list(gene_input)
        assert gene_list == ["smn1"]

        app.genes_to_call = ["smn1"]
        gene_input = "ncf1"
        gene_list = app.get_gene_list(gene_input)
        assert gene_list == ["ncf1"]

        app.genes_to_call = []
        gene_input = None
        gene_list = app.get_gene_list(gene_input)
        assert gene_list == ["smn1", "ncf1"]

        app.genes_to_call = None
        gene_input = None
        gene_list = app.get_gene_list(gene_input)
        assert gene_list == ["smn1", "ncf1"]

        app.genes_to_call = None
        gene_input = "smn1,nfc1"
        try:
            gene_list = app.get_gene_list(gene_input)
        finally:
            assert gene_list == ["smn1"]

        app.genes_to_call = None
        gene_input = "smn,nfc1"
        try:
            gene_list = app.get_gene_list(gene_input)
        finally:
            assert gene_list == []
