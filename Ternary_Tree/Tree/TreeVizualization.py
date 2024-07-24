from __future__ import annotations

from .Nodes import NodeContacts, QubitNum, BranchNum
from .TernaryTree import TernaryTree

from pyvis.network import Network


def mapping(tt: TernaryTree):
    xy = {tt.root: [(0, 0), "root"]}

    def count_qubit_childs(NodeC: NodeContacts):
        return len([child for child in NodeC if not child.is_last])

    def down(node: QubitNum = tt.root, borders: tuple(float) = (-1, 1), deep=1):
        nonlocal xy

        def get_pos(_borders, node, i=3):
            nonlocal xy, deep
            dist = 0.03
            childs = tt[node].childs
            match i:
                case 0:
                    for index, child in enumerate(childs):
                        if isinstance(child, BranchNum):
                            xy[child] = [(dist * (index - 1) + (borders[1] + borders[0]) / 2, deep - 0.5 * 0), node]
                    pass

                case 1:
                    m = (borders[1] + borders[0]) / 2
                    if not childs[0].is_last:
                        xy[childs[0]] = [(m, deep), node]
                        if isinstance(childs[1], BranchNum):
                            xy[childs[1]] = [(m + dist, deep - 0.5 * 0), node]
                        if isinstance(childs[2], BranchNum):
                            xy[childs[2]] = [(m + dist * 2, deep - 0.5 * 0), node]
                    if not childs[1].is_last:
                        xy[childs[1]] = [(m, deep), node]
                        if isinstance(childs[0], BranchNum):
                            xy[childs[0]] = [(m - dist, deep - 0.5 * 0), node]
                        if isinstance(childs[2], BranchNum):
                            xy[childs[2]] = [(m + dist, deep - 0.5 * 0), node]
                    if not childs[2].is_last:
                        xy[childs[2]] = [(m, deep), node]
                        if isinstance(childs[0], BranchNum):
                            xy[childs[0]] = [(m - dist * 2, deep - 0.5 * 0), node]
                        if isinstance(childs[1], BranchNum):
                            xy[childs[1]] = [(m - dist, deep - 0.5 * 0), node]
                    pass

                case 2:
                    l = (borders[1] - borders[0]) / 2
                    j = 0
                    for child in childs:
                        if not child.is_last:
                            _borders = (borders[0] + j * l, borders[0] + (j + 1) * l)
                            xy[child] = [((_borders[1] + _borders[0]) / 2, deep), node]
                            j += 1
                    if isinstance(childs[0], BranchNum):
                        xy[childs[0]] = [((borders[1] + 3 * borders[0]) / 4 - dist, deep - 0.5 * 0), node]
                    if isinstance(childs[1], BranchNum):
                        xy[childs[1]] = [((borders[1] + borders[0]) / 2, deep - 0.5 * 0), node]
                    if isinstance(childs[2], BranchNum):
                        xy[childs[2]] = [((3 * borders[1] + borders[0]) / 4 - dist, deep - 0.5 * 0), node]
                    pass

                case 3:
                    l = (borders[1] - borders[0]) / 3
                    for index, child in enumerate(childs):
                        _borders = (borders[0] + index * l, borders[0] + (index + 1) * l)
                        xy[child] = [((_borders[1] + _borders[0]) / 2, deep), node]

        j = 0
        i = count_qubit_childs(tt[node])
        get_pos(borders, node, i)
        for index, child in enumerate(tt[node]):
            if not child.is_last:
                l = (borders[1] - borders[0]) / i
                _borders = (borders[0] + j * l, borders[0] + (j + 1) * l)
                down(child,
                     borders=_borders,
                     deep=deep + 1)
                j += 1

    down(tt.root)
    return xy


def draw(tt: TernaryTree, k=1):
    xy = mapping(tt)

    from pyvis.network import Network
    import pandas as pd

    got_net = Network(height="1000px", width="100%", bgcolor="#222222", font_color="white")

    # set the physics layout of the network
    # got_net.barnes_hut()

    for node in xy:
        if not node.is_last:
            got_net.add_node(str(node), size=10, fixed=True, title=str(node),
                             x=xy[node][0][0] * 2000 * k, y=xy[node][0][1] * 100)
        else:
            got_net.add_node(str(node), size=10, fixed=True, title=str(node),
                             x=xy[node][0][0] * 2000 * k, y=xy[node][0][1] * 100, color="red")

    for node in xy:
        if isinstance(node, QubitNum):
            if not tt.is_root(node):
                got_net.add_edge(str(node), str(xy[node][1]))
        else:
            got_net.add_edge(str(node), str(xy[node][1]))

    got_net.show("tree.html", notebook=False)
