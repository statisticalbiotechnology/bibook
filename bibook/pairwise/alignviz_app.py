"""A small Streamlit app for exploring pairwise alignment matrices.

Run it with::

    streamlit run bibook/pairwise/alignviz_app.py

The sequences, the scoring function and the alignment mode are widgets; the
dynamic programming matrix is redrawn as you change them, either finished or
stopped at any point during the filling in. The drawing is done by
``alignviz.py``, the same module the notebook uses, and every figure can be
downloaded as an svg.
"""

import base64

import streamlit as st

from alignviz import Alignment, Scoring, Style


def show_svg(svg: str) -> None:
    """Put an svg string on the page. ``st.image`` cannot take one."""
    data = base64.b64encode(svg.encode()).decode()
    st.markdown(
        f'<img src="data:image/svg+xml;base64,{data}" style="max-width:100%">',
        unsafe_allow_html=True,
    )

MODES = {
    "Global (Needleman-Wunsch)": "global",
    "Local (Smith-Waterman)": "local",
    "Semi-global": "semi",
}

st.set_page_config(page_title="Alignment matrices", layout="wide")
st.title("Stages of a pairwise alignment")

with st.sidebar:
    st.header("Sequences")
    seqA = st.text_input("Sequence $a$ (rows)", "TGCATTA").strip().upper()
    seqB = st.text_input("Sequence $b$ (columns)", "GCATTAC").strip().upper()

    st.header("Scoring")
    match = st.number_input("Match", value=3.0, step=1.0)
    mismatch = st.number_input("Mismatch", value=-1.0, step=1.0)
    gap = st.number_input("Gap penalty", value=-2.0, step=1.0)

    st.header("Alignment")
    mode = MODES[st.radio("Type", list(MODES), label_visibility="collapsed")]

    st.header("Drawing")
    cell = st.slider("Cell size", 30, 110, 62)
    show_trace = st.checkbox("Pointers (red)", True)
    show_path = st.checkbox("Traceback (blue)", True)
    show_sources = st.checkbox("Shade the three source cells", True)

if not seqA or not seqB:
    st.info("Give both sequences to get started.")
    st.stop()

aln = Alignment(
    seqA, seqB, scoring=Scoring(match, mismatch, gap), mode=mode
)
style = Style(cell=cell, value_size=cell * 0.35, letter_size=cell * 0.39)

n_steps = len(aln.fill_cells)
step = st.slider(
    "Cells filled in",
    0,
    n_steps,
    n_steps,
    help="Drag back to watch the matrix being filled in, row by row.",
)
finished = step == n_steps

current = None if finished else aln.fill_cells[step]
caption = None
if current is not None:
    i, j = current
    caption = (
        f"$S_{{{i},{j}}}$, comparing {seqA[i - 1]} and {seqB[j - 1]}:  "
        + aln.steps[current].expression()
    )

figure = aln.figure(
    upto=aln.n_init + step,
    trace=show_trace,
    path=aln.traceback() if (finished and show_path) else None,
    current=current,
    sources=show_sources,
    style=style,
)

left, right = st.columns([3, 1])

with left:
    if caption:
        st.markdown(caption)
    show_svg(figure.svg)

with right:
    st.metric("Score", f"{aln.score:g}")
    outA, outB = aln.aligned()
    st.markdown("**Optimal alignment**")
    bars = "".join("|" if x == y and x != "-" else " " for x, y in zip(outA, outB))
    st.code(f"{outA}\n{bars}\n{outB}", language=None)
    if any(len(v) > 1 for v in aln.trace.values()):
        st.caption(
            "Cells with more than one incoming arrow have several optimal "
            "predecessors, so more than one optimal alignment exists."
        )
    st.download_button(
        "Download svg",
        figure.svg,
        file_name=f"{mode}_{seqA}_{seqB}.svg",
        mime="image/svg+xml",
    )
