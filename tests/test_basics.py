"""Basic tests to check that the main examples work."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dna_features_viewer import (BiopythonTranslator, GraphicFeature,
                                 GraphicRecord, CircularGraphicRecord)
from bokeh.resources import CDN
from bokeh.embed import file_html
from Bio import SeqIO
import numpy as np

example_genbank = os.path.join('tests', 'data', "example_sequence.gb")

def test_by_hand(tmpdir):
    """Test building a GraphicRecord "by hand" """
    features = [
        GraphicFeature(start=5, end=20, strand=+1, color="#ffd700",
                       label="Small feature"),
        GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                       label="Gene 1 with a very long name"),
        GraphicFeature(start=400, end=700, strand=-1, color="#cffccc",
                       label="Gene 2"),
        GraphicFeature(start=600, end=900, strand=+1, color="#ccccff",
                       label="Gene 3"),
    ]

    # PLOT AND EXPORT A LINEAR VIEW OF THE CONSTRUCT
    record = GraphicRecord(sequence_length=1000, features=features)
    record.plot(figure_width=5, with_ruler=False) # lazy, just for coverage
    ax, _ = record.plot(figure_width=5)
    target_file = os.path.join(str(tmpdir), "by_hand.png")
    ax.figure.savefig(target_file)

    # PLOT AND EXPORT A CIRCULAR VIEW OF THE CONSTRUCT
    circular_rec = CircularGraphicRecord(sequence_length=1000,
                                         features=features)
    ax2, _ = circular_rec.plot(figure_width=4)
    ax2.figure.tight_layout()
    target_file = os.path.join(str(tmpdir), "by_hand_circular.png")
    ax2.figure.savefig(target_file, bbox_inches="tight")


def test_from_genbank(tmpdir):
    graphic_record = BiopythonTranslator().translate_record(example_genbank)
    ax, _ = graphic_record.plot(figure_width=10)
    ax.figure.tight_layout()
    target_file = os.path.join(str(tmpdir), "from_genbank.png")
    ax.figure.savefig(target_file)

def test_from_genbank_to_circular(tmpdir):
    translator = BiopythonTranslator()
    graphic_record = translator.translate_record(
        example_genbank, record_class=CircularGraphicRecord)
    ax, _ = graphic_record.plot(figure_width=7)
    ax.figure.tight_layout()
    target_file = os.path.join(str(tmpdir), "from_genbank.png")
    ax.figure.savefig(target_file)


def test_plot_with_gc_content(tmpdir):

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 4), sharex=True)

    # Parse the genbank file, plot annotations
    record = SeqIO.read(example_genbank, "genbank")
    graphic_record = BiopythonTranslator().translate_record(record)
    ax, levels = graphic_record.plot()
    graphic_record.plot(ax=ax1, with_ruler=False)

    # Plot the local GC content
    def plot_local_gc_content(record, window_size, ax):
        def gc_content(seq):
            return 100.0*len([c for c in seq if c in "GC"]) / len(seq)
        yy = [gc_content(record.seq[i:i+window_size])
              for i in range(len(record.seq)-window_size)]
        xx = np.arange(len(record.seq)-window_size)+25
        ax.fill_between(xx, yy, alpha=0.3)
        ax.set_ylabel("GC(%)")

    plot_local_gc_content(record, window_size=50, ax=ax2)

    # Resize the figure to the right height
    target_file = os.path.join(str(tmpdir), "with_plot.png")
    fig.tight_layout()
    fig.savefig(target_file)


def test_plot_with_bokeh(tmpdir):
    gb_record = SeqIO.read(example_genbank, "genbank")
    record = BiopythonTranslator().translate_record(record=gb_record)
    plot = record.plot_with_bokeh(figure_width=8)
    target_file = os.path.join(str(tmpdir), "plot_with_bokeh.html")
    with open(target_file, "w+") as f:
        f.write(file_html(plot, CDN, "Example Sequence"))
    with open(target_file, 'r') as f:
        assert len(f.read()) > 5000

def test_plot_with_bokeh_no_labels(tmpdir):
    """Bokeh has a problem with empty lists of labels."""
    gb_record = SeqIO.read(example_genbank, "genbank")
    record = BiopythonTranslator().translate_record(record=gb_record)
    for feature in record.features:
        feature.label = None
    plot = record.plot_with_bokeh(figure_width=8)
    target_file = os.path.join(str(tmpdir), "plot_with_bokeh.html")
    with open(target_file, "w+") as f:
        f.write(file_html(plot, CDN, "Example Sequence"))
    with open(target_file, 'r') as f:
        assert len(f.read()) > 5000

def test_split_overflowing_features():
    features = [
        GraphicFeature(start=10, end=20, strand=+1, label="a"),
        GraphicFeature(start=40, end=55, strand=+1, label="b"),
        GraphicFeature(start=-20, end=2, strand=+1, label="c")
    ]

    # PLOT AND EXPORT A LINEAR VIEW OF THE CONSTRUCT
    record = GraphicRecord(sequence_length=50, features=features)
    record.split_overflowing_features_circularly()
    new_features_locations_and_labels = sorted([
        (f.start, f.end, f.label)
        for f in record.features
    ])
    assert new_features_locations_and_labels == [
        (0, 2, 'c'),
        (0, 5, 'b'),
        (10, 20, 'a'),
        (30, 49, 'c'),
        (40, 49, 'b')
    ]

def test_cropping():
    features=[
        GraphicFeature(start=5, end=20, strand=+1, color="#ffd700",
                       label="Small feature"),
        GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                       label="Gene 1 with a very long name"),
        GraphicFeature(start=400, end=700, strand=-1, color="#cffccc",
                       label="Gene 2"),
        GraphicFeature(start=600, end=900, strand=+1, color="#ccccff",
                       label="Gene 3"),
    ]


    # PLOT AND EXPORT A LINEAR VIEW OF THE CONSTRUCT
    record = GraphicRecord(sequence_length=1000, features=features)
    cropped_record = record.crop((425, 650))
    assert len(cropped_record.features) == 3


def test_to_biopython_record():
    record = GraphicRecord(sequence_length=50, features=[
        GraphicFeature(start=5, end=20, strand=+1, label="a"),
        GraphicFeature(start=20, end=500, strand=+1, label="b"),
        GraphicFeature(start=400, end=700, strand=-1, label="c")
    ])
    biopython_record = record.to_biopython_record(sequence=50*"A")
    features = sorted([
        (f.location.start, f.location.end, f.qualifiers['label'])
        for f in biopython_record.features
    ])
    assert features == [(5, 20, 'a'), (20, 500, 'b'), (400, 700, 'c')]

def test_sequence_and_translation_plotting():
    from dna_features_viewer import (GraphicFeature, GraphicRecord,
                                 CircularGraphicRecord)
    features = [
        GraphicFeature(start=5, end=10, strand=+1, color="#ffd700",
                       label="bbS-1"),
        GraphicFeature(start=8, end=15, strand=+1, color="#ffcccc",
                       label="CrC")
    ]

    record = GraphicRecord(sequence=7*"ATGC", features=features)
    ax, _ = record.plot(figure_width=5)
    record.plot_sequence(ax)
    record.plot_translation(ax, (8, 23), fontdict={'weight': 'bold'})
def test_translate_feature():
    from dna_features_viewer import BiopythonTranslator
    from Bio.SeqFeature import (SeqFeature, ExactPosition, FeatureLocation)

    feat_start = 2
    feat_end = 8

    loc = FeatureLocation(
            ExactPosition(feat_start),
            ExactPosition(feat_end),
    )
    feat = SeqFeature(
            location = loc,
            type = 'gene',
    )

    g_feature = BiopythonTranslator().translate_feature(feat)

    assert (g_feature.start, g_feature.end) == (feat_start, feat_end)
def test_translate_compoundfeature():
    from dna_features_viewer import BiopythonTranslator
    from Bio.SeqFeature import (SeqFeature, ExactPosition, FeatureLocation,
            CompoundLocation)

    feat_start = 14
    feat_end = 6
    seq_len = 20

    loc_spanning_origin = CompoundLocation(
            parts = [
                    FeatureLocation(
                            ExactPosition(feat_start),
                            ExactPosition(seq_len),
                    ),
                    FeatureLocation(
                            ExactPosition(0),
                            ExactPosition(feat_end),
                    ),
            ],
            operator = 'order',
    )
    feat_spanning_origin = SeqFeature(
            location = loc_spanning_origin,
            type = 'gene',
    )

    g_feature = BiopythonTranslator().translate_feature(feat_spanning_origin)

    assert (g_feature.start, g_feature.end) == (feat_start, feat_end)
