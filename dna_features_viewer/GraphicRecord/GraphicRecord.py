import textwrap

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Alphabet import DNAAlphabet

from .matplotlib_plots import MatplotlibPlottableMixin
from .bokeh_plots import BokehPlottableMixin


class GraphicRecord(MatplotlibPlottableMixin, BokehPlottableMixin):
    """Set of Genetic Features of a same DNA sequence, to be plotted together.

    Parameters
    ----------

    sequence_length
      Length of the DNA sequence, in number of nucleotides

    features
      list of GraphicalFeature objects.

    feature_level_height
      Width in inches of one "level" for feature arrows.

    annotation_height
      Width in inches of one "level" for feature annotations.

    first_index
      Indicates the first index to plot in case the sequence is actually a
      subsequence of a larger one. For instance, if the Graphic record
      represents the segment (400, 420) of a sequence, we will have
      ``first_index=400`` and ``sequence_length=20``.

    plots_indexing
      Indicates which standard to use to show nucleotide indices in the plots.
      If 'biopython', the standard python indexing is used (starting at 0).
      If 'genbank', the indexing follows the Genbank standard (starting at 1).
    
    labels_spacing
      Number of pixels that will "pad" every labels to force some horizontal
      space between two labels or between a label and the borders of a feature. 

    Attributes
    ----------

      default_font_family
        Default font to use for a feature that doesn't declare a font.
      
      default_ruler_color
        Default ruler color to use when no color is given at plot() time.
    
      default_box_color
        Default box color for non-inline annotations. If set to None, no
        boxes will be drawn unless the features declare a box_color.
        If "auto", a color (clearer version of the feature's color) will be
        computed, for all features also declaring their box_color as "auto".
    """
    
    default_font_family = None
    default_ruler_color = 'grey'
    default_box_color = 'auto'
    
    def __init__(
        self,
        sequence_length=None,
        sequence=None,
        features=(),
        feature_level_height=1,
        first_index=0,
        plots_indexing="biopython",
        labels_spacing=8,
    ):
        if sequence_length is None:
            sequence_length = len(sequence)
        self.features = features
        self.sequence_length = sequence_length
        self.feature_level_height = feature_level_height
        self.sequence = sequence
        self.first_index = first_index
        self.plots_indexing = plots_indexing
        self.labels_spacing = labels_spacing

    @property
    def span(self):
        """Return the display span (start, end) accounting for first_index."""
        return self.first_index, self.first_index + self.sequence_length

    def to_biopython_record(self, sequence):
        """
        Example
        -------
        from Bio import SeqIO
        gr_record = GraphicRecord(features=features, sequence_length=len(seq),
                                  sequence=seq)
        bio_record = gr_record.to_biopython_record()
        with open("example.gb", "w+") as f:
            SeqIO.write(record, f, "genbank")
        """
        features = [
            SeqFeature(
                FeatureLocation(f.start, f.end, f.strand),
                type=f.feature_type,
                qualifiers={"label": f.label},
            )
            for f in self.features
        ]
        if not isinstance(sequence, Seq):
            sequence = Seq(sequence, alphabet=DNAAlphabet())
        return SeqRecord(seq=sequence, features=features)

    def crop(self, window):
        s, e = window
        if (s < 0) or (e >= self.sequence_length):
            raise ValueError("out-of-bound cropping")
        new_features = []
        for f in self.features:
            cropped_feature = f.crop(window)
            if cropped_feature is not None:  # = has ovelap with the window
                new_features.append(cropped_feature)

        return GraphicRecord(
            sequence=self.sequence[s:e] if self.sequence is not None else None,
            sequence_length=e - s,
            features=new_features,
            feature_level_height=self.feature_level_height,
            first_index=self.first_index + s,
        )

    def determine_annotation_height(self, levels):
        """By default the ideal annotation level height is the same as the
        feature_level_height."""
        # Improve me: ideally, annotation width would be linked to the height
        # of one line of text, so dependent on font size and ax height/span.
        return self.feature_level_height

    def coordinates_in_plot(self, x, level):
        """Convert a sequence position and height level into a (x, y) position.
        """
        return (x, level * self.feature_level_height)

    def split_overflowing_features_circularly(self):
        """Split the features that overflow over the edge for circular
        constructs (inplace)."""
        new_features = []
        for f in self.features:
            if f.start < 0 < f.end:
                f1, f2 = f.split_in_two(-1)
                f1.start, f1.end = (
                    f1.start + self.sequence_length,
                    f1.end + self.sequence_length,
                )
                new_features += [f1, f2]
            elif f.start < self.sequence_length < f.end:
                f1, f2 = f.split_in_two(self.sequence_length - 1)
                f2.start, f2.end = (
                    f2.start - self.sequence_length,
                    f2.end - self.sequence_length,
                )
                new_features += [f1, f2]
            else:
                new_features.append(f)
        self.features = new_features

    def _format_label(self, label, max_label_length=50, max_line_length=40):
        if len(label) > max_label_length:
            label = label[: max_label_length - 1] + "…"
        label = "\n".join(textwrap.wrap(label, max_line_length))
        return label
    
    def compute_padding(self, ax):
        ax_width = ax.get_window_extent(ax.figure.canvas.get_renderer()).width
        xmin, xmax = ax.get_xlim()
        return self.labels_spacing * (xmax - xmin) / (1.0 * ax_width)
