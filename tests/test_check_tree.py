"""Tests for workflow/check_tree.py"""

import os
import sys
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))

try:
    from ete3 import Tree
    import check_tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False

pytestmark = pytest.mark.skipif(not ETE3_AVAILABLE, reason="ete3 not installed")


# ---------------------------------------------------------------------------
# check_strict_binary
# ---------------------------------------------------------------------------

def test_check_strict_binary_valid():
    # A simple 3-leaf balanced binary tree: ((A,B),C)
    t = Tree("((A,B),C);")
    check_tree.check_strict_binary(t)  # should not raise


def test_check_strict_binary_larger_valid():
    t = Tree("(((A,B),(C,D)),((E,F),(G,H)));")
    check_tree.check_strict_binary(t)  # should not raise


def test_check_strict_binary_trifurcation_raises():
    # Root has 3 children → not strictly binary
    t = Tree("(A,B,C);")
    with pytest.raises(ValueError, match="not strictly binary"):
        check_tree.check_strict_binary(t)


def test_check_strict_binary_internal_polytomy_raises():
    # Internal node with 3 children
    t = Tree("((A,B,C),D);")
    with pytest.raises(ValueError, match="not strictly binary"):
        check_tree.check_strict_binary(t)


def test_check_strict_binary_single_leaf_passes():
    # A tree with a single leaf has no internal nodes → trivially binary
    t = Tree("A;")
    check_tree.check_strict_binary(t)


# ---------------------------------------------------------------------------
# collapse_unary
# ---------------------------------------------------------------------------

def test_collapse_unary_removes_single_child_internal_node():
    """
    Build a tree where an internal node has only one child and verify
    that after collapse_unary, all remaining internal nodes have ≥ 2 children.
    """
    # Manually construct a unary node: root -> internal(only_child=A), B
    t = Tree()
    root = t
    unary = root.add_child()
    leaf_a = unary.add_child(name="A")
    leaf_b = root.add_child(name="B")

    check_tree.collapse_unary(t)

    for node in t.traverse():
        if not node.is_leaf():
            assert len(node.children) != 1, (
                f"Node '{node.name}' still has exactly 1 child after collapse_unary"
            )


def test_collapse_unary_no_op_on_binary_tree():
    t = Tree("((A,B),(C,D));")
    original_leaves = {leaf.name for leaf in t.iter_leaves()}

    check_tree.collapse_unary(t)

    # Leaves must be preserved
    assert {leaf.name for leaf in t.iter_leaves()} == original_leaves
    # All internal nodes still have 2 children
    for node in t.traverse():
        if not node.is_leaf():
            assert len(node.children) == 2
