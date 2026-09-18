"""SVG illustrations of dynamic programming matrices for pairwise alignment.

This module is an authoring tool for the book. It replaces the TikZ/LaTeX
route (see ``nw_code_latex.ipynb``) with plain SVG, so that figures can be
produced, inspected and tweaked directly inside a notebook without a LaTeX
installation.

The typical use is::

    from alignviz import Alignment

    aln = Alignment("GAC", "ACG", match=1, mismatch=-1, gap=-1)
    aln.figure(upto=aln.n_init)          # the initialised borders
    aln.figure()                         # the completed matrix
    aln.figure(path=aln.traceback())     # with the optimal path

Every figure renders itself in Jupyter, and can be written to disk with
``figure.save("img/nw_short_1.svg")``.

The three classical variants are supported through ``mode``:

===========  =====================  =========================  ====================
``mode``     initialisation         recursion                   traceback starts at
===========  =====================  =========================  ====================
``global``   cumulative gap costs   max of the three moves      bottom right cell
``local``    zeros                  max of the three moves, 0   the largest cell
``semi``     zeros                  max of the three moves      max of last row/col
===========  =====================  =========================  ====================
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = ["Scoring", "Alignment", "Figure", "Style", "player"]

Cell = Tuple[int, int]
Move = Tuple[int, int]

DIAG: Move = (-1, -1)
UP: Move = (-1, 0)
LEFT: Move = (0, -1)

MOVE_NAMES = {DIAG: "diagonal", UP: "up", LEFT: "left"}

# Which corner of a cell each of the three competing values is written in,
# as (right, down) multipliers away from the middle of the cell. A candidate
# is placed towards the cell it comes from: the diagonal predecessor lies up
# and to the left, so it goes in the top left corner. The other two would
# both want the corner they point at, so the one from above takes the top
# right corner and the one from the left the bottom left one.
CANDIDATE_CORNERS: Dict[Move, Move] = {DIAG: (-1, -1), UP: (1, -1), LEFT: (-1, 1)}


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------


@dataclass
class Scoring:
    """A simple match/mismatch/gap scoring function, ``d(x,y)``."""

    match: float = 3.0
    mismatch: float = -1.0
    gap: float = -2.0

    def __call__(self, a: str, b: str) -> float:
        if a == "-" or b == "-":
            return self.gap
        return self.match if a == b else self.mismatch


# --------------------------------------------------------------------------
# Drawing style
# --------------------------------------------------------------------------


@dataclass
class Style:
    """Everything that decides how a matrix looks."""

    cell: float = 62.0  # side of a matrix cell, in px
    pad: float = 6.0  # margin around the whole drawing
    font_family: str = "Latin Modern Roman, Computer Modern, Georgia, serif"
    value_size: float = 22.0  # font size of the matrix values
    letter_size: float = 24.0  # font size of the sequence letters
    grid_color: str = "#000000"
    grid_width: float = 1.2
    text_color: str = "#000000"
    trace_color: str = "#e00000"  # the red "how did I get here" arrows
    trace_width: float = 2.6
    trace_dash: str = "5,4"
    path_color: str = "#0033cc"  # the blue traceback arrows
    path_width: float = 3.2
    current_fill: str = "#fff3bf"  # cell being computed
    source_fill: str = "#dbe9ff"  # cells it is computed from
    candidate_size: Optional[float] = None  # the competing values, in the
    candidate_color: str = "#666666"  # corners of the cell being computed
    mark_fill: str = "#e8f7e8"  # user supplied highlights
    caption_size: float = 20.0
    caption_color: str = "#222222"


# --------------------------------------------------------------------------
# The dynamic programming matrix
# --------------------------------------------------------------------------


@dataclass
class Step:
    """The bookkeeping of a single application of the recursion."""

    i: int
    j: int
    candidates: List[Tuple[Move, float, float]]  # move, previous score, increment
    value: float
    moves: List[Move]

    def expression(self, fmt=None) -> str:
        """``max(0+-1,-1+-1,-1+-1)=-1`` - the caption used in the book."""
        fmt = fmt or _fmt
        terms = ",".join(f"{fmt(prev)}+{fmt(inc)}" for _, prev, inc in self.candidates)
        return f"max({terms})={fmt(self.value)}"


def _fmt(x: float) -> str:
    """Integers without decimals, everything else with one decimal."""
    if x == int(x):
        return str(int(x))
    return f"{x:.1f}"


class Alignment:
    """Fill in a dynamic programming matrix and draw it at any stage."""

    def __init__(
        self,
        seqA: str,
        seqB: str,
        scoring: Optional[Scoring] = None,
        match: Optional[float] = None,
        mismatch: Optional[float] = None,
        gap: Optional[float] = None,
        mode: str = "global",
    ):
        if mode not in ("global", "local", "semi"):
            raise ValueError("mode must be 'global', 'local' or 'semi'")
        if scoring is None:
            scoring = Scoring()
        if match is not None:
            scoring.match = float(match)
        if mismatch is not None:
            scoring.mismatch = float(mismatch)
        if gap is not None:
            scoring.gap = float(gap)

        self.seqA, self.seqB = seqA, seqB
        self.scoring, self.mode = scoring, mode
        self.m, self.n = len(seqA) + 1, len(seqB) + 1

        self.S: List[List[float]] = [[0.0] * self.n for _ in range(self.m)]
        self.trace: Dict[Cell, List[Move]] = {}
        self.steps: Dict[Cell, Step] = {}

        self._initiate()
        self._fill()

        # The order in which the cells are filled: first the borders (which
        # count as one single step, they are all set by the initialisation),
        # then the interior of the matrix, row by row.
        self.init_cells: List[Cell] = (
            [(0, 0)]
            + [(i, 0) for i in range(1, self.m)]
            + [(0, j) for j in range(1, self.n)]
        )
        self.fill_cells: List[Cell] = [
            (i, j) for i in range(1, self.m) for j in range(1, self.n)
        ]
        self.order: List[Cell] = self.init_cells + self.fill_cells
        self.n_init = len(self.init_cells)

    # -- the algorithm ------------------------------------------------------

    def _initiate(self) -> None:
        g = self.scoring.gap
        self.trace[(0, 0)] = []
        for i in range(1, self.m):
            self.S[i][0] = i * g if self.mode == "global" else 0.0
            self.trace[(i, 0)] = [UP] if self.mode == "global" else []
        for j in range(1, self.n):
            self.S[0][j] = j * g if self.mode == "global" else 0.0
            self.trace[(0, j)] = [LEFT] if self.mode == "global" else []

    def _fill(self) -> None:
        d = self.scoring
        for i in range(1, self.m):
            for j in range(1, self.n):
                a, b = self.seqA[i - 1], self.seqB[j - 1]
                candidates = [
                    (DIAG, self.S[i - 1][j - 1], d(a, b)),
                    (UP, self.S[i - 1][j], d(a, "-")),
                    (LEFT, self.S[i][j - 1], d("-", b)),
                ]
                best = max(prev + inc for _, prev, inc in candidates)
                if self.mode == "local":
                    best = max(best, 0.0)
                self.S[i][j] = best
                moves = [mv for mv, prev, inc in candidates if prev + inc == best]
                self.trace[(i, j)] = moves
                self.steps[(i, j)] = Step(i, j, candidates, best, moves)

    # -- results ------------------------------------------------------------

    @property
    def score(self) -> float:
        """The score of the optimal alignment, for the chosen mode."""
        return self.S[self.start()[0]][self.start()[1]]

    def start(self) -> Cell:
        """The cell the traceback starts from; the first of ``starts()``."""
        return self.starts()[0]

    def starts(self) -> List[Cell]:
        """Every cell an optimal traceback may start from.

        The bottom right corner for a global alignment, every cell holding
        the best score for a local one, and every cell of the last row or
        column that does for a semi-global one. There can be more than one,
        which is one of the two ways several optimal alignments arise.
        """
        if self.mode == "global":
            return [(self.m - 1, self.n - 1)]
        if self.mode == "local":
            cells = [(i, j) for i in range(self.m) for j in range(self.n)]
        else:
            cells = [(self.m - 1, j) for j in range(self.n)]
            cells += [(i, self.n - 1) for i in range(self.m)]
            cells = list(dict.fromkeys(cells))  # the corner is in both lists
        best = max(self.S[i][j] for i, j in cells)
        if self.mode == "local" and best <= 0:
            # Nothing beats the empty alignment, so do not call every zero
            # in the matrix a starting point.
            return [max(cells, key=lambda c: self.S[c[0]][c[1]])]
        return [c for c in cells if self.S[c[0]][c[1]] == best]

    def _stop(self, cell: Cell) -> bool:
        i, j = cell
        if self.mode == "global":
            return i == 0 and j == 0
        if self.mode == "local":
            return self.S[i][j] <= 0
        return i == 0 or j == 0

    def traceback(self, start: Optional[Cell] = None) -> List[Cell]:
        """One optimal path, as a list of cells from the start cell backwards.

        Where several moves are optimal the first one (diagonal before up
        before left) is taken, which is the same convention as in the book.
        """
        i, j = start or self.start()
        path = [(i, j)]
        while not self._stop((i, j)):
            moves = self.trace.get((i, j), [])
            if not moves:
                break
            di, dj = moves[0]
            i, j = i + di, j + dj
            path.append((i, j))
        return path

    def tracebacks(
        self, start: Optional[Cell] = None, limit: Optional[int] = None
    ) -> List[List[Cell]]:
        """Every optimal path, not just the first one.

        Wherever a cell has several optimal predecessors the traceback
        branches, and each branch spells out an alignment of its own. Pass
        ``limit`` to stop after that many paths; ``count_tracebacks`` says
        how many there are in total, which can be a very large number.
        """
        out: List[List[Cell]] = []

        def walk(path: List[Cell]) -> None:
            if limit is not None and len(out) >= limit:
                return
            cell = path[-1]
            moves = self.trace.get(cell, [])
            if self._stop(cell) or not moves:
                out.append(path)
                return
            for di, dj in moves:
                walk(path + [(cell[0] + di, cell[1] + dj)])
                if limit is not None and len(out) >= limit:
                    return

        for cell in [start] if start else self.starts():
            walk([cell])
        return out

    def count_tracebacks(self, start: Optional[Cell] = None) -> int:
        """How many optimal paths there are, without enumerating them."""
        memo: Dict[Cell, int] = {}

        def count(cell: Cell) -> int:
            moves = self.trace.get(cell, [])
            if self._stop(cell) or not moves:
                return 1
            if cell not in memo:
                memo[cell] = sum(
                    count((cell[0] + di, cell[1] + dj)) for di, dj in moves
                )
            return memo[cell]

        return sum(count(c) for c in ([start] if start else self.starts()))

    def aligned(self, path: Optional[List[Cell]] = None) -> Tuple[str, str]:
        """The two aligned sequences that a traceback path spells out."""
        path = path or self.traceback()
        outA, outB = "", ""
        for (i, j), (pi, pj) in zip(path, path[1:]):
            outA = (self.seqA[pi] if i - pi else "-") + outA
            outB = (self.seqB[pj] if j - pj else "-") + outB
        return outA, outB

    def alignment_score(self, path: Optional[List[Cell]] = None) -> float:
        i, j = (path or self.traceback())[0]
        return self.S[i][j]

    # -- drawing ------------------------------------------------------------

    def figure(
        self,
        upto: Optional[int] = None,
        trace: bool = True,
        path: Optional[Sequence[Cell]] = None,
        paths: Optional[Sequence[Sequence[Cell]]] = None,
        visiting: Optional[Iterable[Cell]] = None,
        current: Optional[Cell] = None,
        sources: bool = False,
        candidates: bool = False,
        highlight: Optional[Iterable[Cell]] = None,
        caption: Optional[str] = None,
        style: Optional[Style] = None,
        fmt=None,
    ) -> "Figure":
        """Draw the matrix as it looks after ``upto`` cells have been filled.

        ``upto`` counts cells in ``self.order``; ``self.n_init`` therefore
        gives the freshly initialised matrix and ``None`` the finished one.
        Set ``current`` to a cell to highlight it, ``sources=True`` to
        also shade the (up to) three cells the recursion reads from, and
        ``candidates=True`` to write the three competing values in the
        corners of the current cell, each in the corner nearest the cell it
        comes from.

        ``path`` draws one traceback and ``paths`` any number of them, with
        the arrows shared by several paths drawn only once. Truncated paths
        draw a traceback that is still under way; ``visiting`` then shades
        the cells it has reached so far.
        """
        style = style or Style()
        fmt = fmt or _fmt
        if upto is None:
            upto = len(self.order)
        shown = set(self.order[:upto])

        marks: Dict[Cell, str] = {}
        for cell in highlight or ():
            marks[cell] = style.mark_fill
        drawn = [list(q) for q in paths] if paths else []
        if path:
            drawn.append(list(path))

        if current is not None and sources:
            i, j = current
            for di, dj in (DIAG, UP, LEFT):
                if 0 <= i + di and 0 <= j + dj:
                    marks[(i + di, j + dj)] = style.source_fill
        for cell in visiting or ():
            marks[cell] = style.current_fill
        if current is not None:
            marks[current] = style.current_fill

        svg = _render(
            self,
            shown=shown,
            trace=trace,
            paths=drawn,
            marks=marks,
            candidate_cell=current if candidates else None,
            caption=caption,
            style=style,
            fmt=fmt,
        )
        return Figure(svg, caption)

    def frames(
        self,
        sources: bool = True,
        candidates: bool = False,
        traceback: bool = True,
        style: Optional[Style] = None,
        fmt=None,
    ) -> List["Figure"]:
        """The whole story: initialisation, one frame per cell, traceback.

        Each frame carries the recursion it just evaluated as its caption,
        which is exactly the text printed under the figures in the book.
        """
        fmt = fmt or _fmt
        out = [
            self.figure(
                upto=self.n_init,
                caption="Initialisation of the borders",
                style=style,
                fmt=fmt,
            )
        ]
        for k, cell in enumerate(self.fill_cells, start=1):
            out.append(
                self.figure(
                    upto=self.n_init + k,
                    current=cell,
                    sources=sources,
                    candidates=candidates,
                    caption=self.steps[cell].expression(fmt),
                    style=style,
                    fmt=fmt,
                )
            )
        if traceback:
            path = self.traceback()
            a, b = self.aligned(path)
            out.append(
                self.figure(
                    path=path,
                    caption=f"Traceback: {a} / {b}, score {fmt(self.alignment_score(path))}",
                    style=style,
                    fmt=fmt,
                )
            )
        return out

    # -- convenience --------------------------------------------------------

    def save_frames(self, pattern: str, **kwargs) -> List[str]:
        """Write every frame to disk. ``pattern`` must contain a ``{}``.

        ``aln.save_frames("img/nw_short_{}.svg")`` reproduces the numbering
        used by the figures in the Needleman-Wunsch chapter.
        """
        paths = []
        for k, fig in enumerate(self.frames(**kwargs), start=1):
            path = pattern.format(k)
            fig.save(path)
            paths.append(path)
        return paths

    def _repr_html_(self) -> str:
        return self.figure()._repr_html_()

    def __repr__(self) -> str:
        return (
            f"Alignment({self.seqA!r}, {self.seqB!r}, mode={self.mode!r}, "
            f"score={_fmt(self.score)})"
        )


# --------------------------------------------------------------------------
# The figure object
# --------------------------------------------------------------------------


class Figure:
    """An SVG drawing that knows how to show itself and how to be saved."""

    def __init__(self, svg: str, caption: Optional[str] = None):
        self.svg = svg
        self.caption = caption

    def _repr_html_(self) -> str:
        return self.svg

    def _repr_svg_(self) -> str:
        return self.svg

    def save(self, path: str) -> str:
        """Write to ``.svg``, or to ``.png``/``.pdf`` if cairosvg is around."""
        directory = os.path.dirname(path)
        if directory:
            os.makedirs(directory, exist_ok=True)
        ext = os.path.splitext(path)[1].lower()
        if ext == ".svg":
            with open(path, "w") as fh:
                fh.write(self.svg)
        elif ext in (".png", ".pdf"):
            try:
                import cairosvg
            except ImportError as exc:  # pragma: no cover
                raise ImportError(
                    f"writing {ext} needs cairosvg (pip install cairosvg); "
                    "saving as .svg always works"
                ) from exc
            writer = cairosvg.svg2png if ext == ".png" else cairosvg.svg2pdf
            writer(bytestring=self.svg.encode(), write_to=path, scale=2)
        else:
            raise ValueError(f"cannot write {ext} files")
        return path

    def __repr__(self) -> str:
        return f"<Figure {len(self.svg)} bytes>"


# --------------------------------------------------------------------------
# SVG generation
# --------------------------------------------------------------------------


def _esc(text: str) -> str:
    return (
        str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    )


def _render(
    aln: Alignment,
    shown: set,
    trace: bool,
    paths: List[List[Cell]],
    marks: Dict[Cell, str],
    candidate_cell: Optional[Cell],
    caption: Optional[str],
    style: Style,
    fmt,
) -> str:
    c = style.cell
    m, n = aln.m, aln.n
    # One cell of margin on the top and on the left carries the sequences.
    x0, y0 = style.pad + c, style.pad + c
    width = 2 * style.pad + c * (n + 1)
    cap_h = (style.caption_size * 1.8) if caption else 0.0
    height = 2 * style.pad + c * (m + 1) + cap_h

    def cx(j: float) -> float:
        return x0 + (j + 0.5) * c

    def cy(i: float) -> float:
        return y0 + (i + 0.5) * c

    out: List[str] = []
    out.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{_num(width)}" '
        f'height="{_num(height)}" viewBox="0 0 {_num(width)} {_num(height)}" '
        f'font-family="{style.font_family}">'
    )
    out.append(f'<rect width="{_num(width)}" height="{_num(height)}" fill="white"/>')

    # Shaded cells, drawn below the grid lines.
    for (i, j), colour in marks.items():
        if 0 <= i < m and 0 <= j < n:
            out.append(
                f'<rect x="{_num(x0 + j * c)}" y="{_num(y0 + i * c)}" '
                f'width="{_num(c)}" height="{_num(c)}" fill="{colour}"/>'
            )

    # The grid.
    g = f'stroke="{style.grid_color}" stroke-width="{_num(style.grid_width)}"'
    for i in range(m + 1):
        y = y0 + i * c
        out.append(f'<line x1="{_num(x0)}" y1="{_num(y)}" x2="{_num(x0 + n * c)}" y2="{_num(y)}" {g}/>')
    for j in range(n + 1):
        x = x0 + j * c
        out.append(f'<line x1="{_num(x)}" y1="{_num(y0)}" x2="{_num(x)}" y2="{_num(y0 + m * c)}" {g}/>')

    # The sequences, along the top and the left hand side.
    for i, letter in enumerate("-" + aln.seqA):
        out.append(_text(cx(-1), cy(i), letter, style.letter_size, style.text_color))
    for j, letter in enumerate("-" + aln.seqB):
        out.append(_text(cx(j), cy(-1), letter, style.letter_size, style.text_color))

    # The values of the cells that have been computed so far.
    for (i, j) in sorted(shown):
        out.append(
            _text(cx(j), cy(i), fmt(aln.S[i][j]), style.value_size, style.text_color)
        )

    # The three competing values, in the corners of the cell being computed.
    step = aln.steps.get(candidate_cell) if candidate_cell is not None else None
    if step is not None:
        size = style.candidate_size or style.value_size * 0.55
        off = 0.30 * c
        for mv, prev, inc in step.candidates:
            sx, sy = CANDIDATE_CORNERS[mv]
            out.append(
                _text(
                    cx(step.j) + sx * off,
                    cy(step.i) + sy * off,
                    fmt(prev + inc),
                    size,
                    style.candidate_color,
                )
            )

    # Red dashed arrows: which cell did the value come from?
    if trace:
        for (i, j) in sorted(shown):
            for mv in aln.trace.get((i, j), []):
                if (i + mv[0], j + mv[1]) in shown:
                    out.append(_trace_arrow(i, j, mv, cx, cy, c, style))

    # Blue arrows: the optimal paths, drawn backwards. Where several paths
    # run through the same cells the arrow is only drawn once.
    seen: set = set()
    for p in paths:
        for (i, j), (pi, pj) in zip(p, p[1:]):
            move = (pi - i, pj - j)
            if ((i, j), move) in seen:
                continue
            seen.add(((i, j), move))
            out.append(_path_arrow(i, j, move, cx, cy, c, style))

    if caption:
        out.append(
            _text(
                width / 2,
                height - style.pad - style.caption_size * 0.55,
                caption,
                style.caption_size,
                style.caption_color,
            )
        )

    out.append("</svg>")
    return "".join(out)


def _num(x: float) -> str:
    return f"{x:.2f}".rstrip("0").rstrip(".")


def _arrowhead(x: float, y: float, ux: float, uy: float,
               colour: str, width: float) -> str:
    """A solid triangular arrowhead at ``(x, y)``, pointing along ``(ux, uy)``.

    Drawn as an explicit polygon rather than as an SVG ``<marker>``. A marker
    is referenced by id, and the book inlines a dozen of these drawings into a
    single html page, where every ``url(#...)`` would resolve to the first
    definition in the document - which may well sit inside a hidden cell, and
    then paints nothing at all.
    """
    length, half = 5.0 * width, 2.5 * width
    px, py = -uy, ux  # perpendicular to the direction of travel
    bx, by = x - ux * length, y - uy * length
    points = " ".join(
        f"{_num(a)},{_num(b)}"
        for a, b in (
            (x, y),
            (bx + px * half, by + py * half),
            (bx - px * half, by - py * half),
        )
    )
    return f'<polygon points="{points}" fill="{colour}"/>'


def _text(x: float, y: float, s: str, size: float, colour: str) -> str:
    return (
        f'<text x="{_num(x)}" y="{_num(y)}" font-size="{_num(size)}" '
        f'fill="{colour}" text-anchor="middle" dominant-baseline="central">'
        f"{_esc(s)}</text>"
    )


def _trace_arrow(i, j, move, cx, cy, c, style: Style) -> str:
    """A short dashed arrow across the border into cell ``(i,j)``."""
    di, dj = move
    # The point on the border between the cell and its predecessor.
    bx = cx(j) + dj * c / 2
    by = cy(i) + di * c / 2
    # Unit vector pointing from the predecessor into (i, j).
    ux, uy = -dj, -di
    norm = math.hypot(ux, uy)
    ux, uy = ux / norm, uy / norm
    reach = 0.23 * c
    ex, ey = bx + ux * reach, by + uy * reach
    return (
        f'<line x1="{_num(bx - ux * reach)}" y1="{_num(by - uy * reach)}" '
        f'x2="{_num(ex)}" y2="{_num(ey)}" '
        f'stroke="{style.trace_color}" stroke-width="{_num(style.trace_width)}" '
        f'stroke-dasharray="{style.trace_dash}"/>'
        + _arrowhead(ex, ey, ux, uy, style.trace_color, style.trace_width)
    )


def _path_arrow(i, j, move, cx, cy, c, style: Style) -> str:
    """A solid arrow from cell ``(i,j)`` back to its predecessor."""
    di, dj = move
    bx = cx(j) + dj * c / 2
    by = cy(i) + di * c / 2
    # Unit vector pointing backwards, from (i, j) towards the predecessor.
    ux, uy = dj, di
    norm = math.hypot(ux, uy)
    ux, uy = ux / norm, uy / norm
    # Offset sideways so that the blue arrow does not sit on top of the red one.
    off = 0.14 * c
    ox, oy = -uy * off, ux * off
    reach = 0.32 * c
    ex, ey = bx + ux * reach + ox, by + uy * reach + oy
    return (
        f'<line x1="{_num(bx - ux * reach + ox)}" y1="{_num(by - uy * reach + oy)}" '
        f'x2="{_num(ex)}" y2="{_num(ey)}" '
        f'stroke="{style.path_color}" stroke-width="{_num(style.path_width)}" '
        f'stroke-linecap="round"/>'
        + _arrowhead(ex, ey, ux, uy, style.path_color, style.path_width)
    )


# --------------------------------------------------------------------------
# A little step-by-step player, for notebooks and for the built book
# --------------------------------------------------------------------------

_PLAYER_TEMPLATE = """
<div class="alignviz" id="{uid}" style="max-width:100%;">
  <div class="alignviz-stage" style="min-height:{stage}px;">{frames}</div>
  <div style="margin-top:.5em;display:flex;gap:.5em;align-items:center;
              font-family:sans-serif;font-size:14px;">
    <button type="button" data-av="first">&#124;&#9664;</button>
    <button type="button" data-av="prev">&#9664;</button>
    <button type="button" data-av="play">&#9654; play</button>
    <button type="button" data-av="next">&#9654;</button>
    <button type="button" data-av="last">&#9654;&#124;</button>
    <input type="range" min="0" max="{last}" value="0" data-av="slider"
           style="flex:1;min-width:80px;">
    <span data-av="label">1/{count}</span>
  </div>
</div>
<script>
(function() {{
  var root = document.getElementById("{uid}");
  var frames = root.querySelectorAll(".alignviz-frame");
  var slider = root.querySelector('[data-av="slider"]');
  var label = root.querySelector('[data-av="label"]');
  var timer = null, at = 0;
  function show(k) {{
    at = Math.max(0, Math.min(frames.length - 1, k));
    for (var i = 0; i < frames.length; i++)
      frames[i].style.display = (i === at) ? "block" : "none";
    slider.value = at;
    label.textContent = (at + 1) + "/" + frames.length;
  }}
  function stop() {{
    if (timer) {{ clearInterval(timer); timer = null; }}
    root.querySelector('[data-av="play"]').innerHTML = "&#9654; play";
  }}
  root.querySelector('[data-av="first"]').onclick = function() {{ stop(); show(0); }};
  root.querySelector('[data-av="prev"]').onclick = function() {{ stop(); show(at - 1); }};
  root.querySelector('[data-av="next"]').onclick = function() {{ stop(); show(at + 1); }};
  root.querySelector('[data-av="last"]').onclick = function() {{ stop(); show(frames.length - 1); }};
  root.querySelector('[data-av="play"]').onclick = function() {{
    if (timer) {{ stop(); return; }}
    this.innerHTML = "&#10074;&#10074; pause";
    timer = setInterval(function() {{
      if (at >= frames.length - 1) {{ stop(); return; }}
      show(at + 1);
    }}, {delay});
  }};
  slider.oninput = function() {{ stop(); show(parseInt(this.value, 10)); }};
  show(0);
}})();
</script>
"""

_uid_counter = [0]


def player(figures: Sequence[Figure], delay: int = 700, height: Optional[int] = None):
    """Wrap a list of figures in a self contained step-through widget.

    No kernel and no javascript library is needed, so the result also works
    in the built html book. Returns an object that renders in Jupyter.
    """
    _uid_counter[0] += 1
    uid = f"alignviz{_uid_counter[0]}"
    frames = "".join(
        f'<div class="alignviz-frame" style="display:none;">{fig.svg}</div>'
        for fig in figures
    )
    if height is None:
        height = _svg_height(figures[0].svg) if figures else 200
    html = _PLAYER_TEMPLATE.format(
        uid=uid,
        frames=frames,
        last=max(0, len(figures) - 1),
        count=len(figures),
        delay=int(delay),
        stage=int(height) + 10,
    )
    return _Html(html)


def _svg_height(svg: str) -> float:
    marker = 'height="'
    start = svg.index(marker, svg.index("<svg")) + len(marker)
    return float(svg[start : svg.index('"', start)])


class _Html:
    def __init__(self, html: str):
        self.html = html

    def _repr_html_(self) -> str:
        return self.html

    def save(self, path: str) -> str:
        with open(path, "w") as fh:
            fh.write(f"<!doctype html><meta charset='utf-8'>{self.html}")
        return path
