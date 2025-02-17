from __future__ import annotations


from abc import abstractmethod, ABC


ROOT = "root"


class NodeNum(ABC):
    def __init__(self, num: int = 0):
        self._num = num
        self.num = self._num

    @property
    @abstractmethod
    def num(self):
        return self._num

    @num.setter
    @abstractmethod
    def num(self, value):
        pass

    @property
    @abstractmethod
    def name(self):
        pass

    @staticmethod
    @abstractmethod
    def _get_num(other):
        pass

    @abstractmethod
    def is_last(self):
        pass

    def __str__(self):
        return str(self._num) + self.name

    def __repr__(self):
        return str(self._num) + self.name

    def __bool__(self):
        return True

    def __eq__(self, other):
        other = self._get_num(other)
        if other is None:
            return False
        else:
            return self.num == other

    def __ne__(self ,other):
        other = self._get_num(other)
        if other is None:
            return False
        else:
            return self.num != other

    def __gt__(self, other):
        other = self._get_num(other)
        if other is None:
            return False
        else:
            return self.num > other

    def __lt__(self, other):
        other = self._get_num(other)
        if other is None:
            return False
        else:
            return self.num < other

    def __hash__(self):
        return hash((self.num, self.name))


class QubitNum(NodeNum):
    """
        Type for numeration of the branches. This is special type for tree leaves.
    """
    def __init__(self, num: int = 0):
        super().__init__(num)

    @property
    def num(self):
        return self._num

    @num.setter
    def num(self, value):
        value = self._get_num(value)
        if value is None:
            raise ValueError("type value shold be int or NodeNum not " + str(type(value)))
        if value >= 0:
            self._num = value
        else:
            self._num = 0

    @property
    def name(self):
        return "q"

    @property
    def is_last(self):
        return False

    @staticmethod
    def _get_num(other):
        if isinstance(other, QubitNum):
            return other.num
        if isinstance(other, int):
            return other
        return None


class BranchNum(NodeNum):
    """
        Type for numeration of the branches. This is special type for tree leaves.
    """
    def __init__(self, num: int = 1):
        super().__init__(num)
        self.sign = +1

    @property
    def num(self):
        return self._num

    @num.setter
    def num(self, value):
        value = self._get_num(value)
        if value > 0:
            self._num = value

    @property
    def name(self):
        return "b"

    @property
    def is_last(self):
        return True

    @staticmethod
    def _get_num(other):
        if isinstance(other, BranchNum):
            return other.num
        if isinstance(other, int):
            return other
        # raise Exception("Attempt to use numeration from" + str(type(other)) + " for BranchNum class")
        return None


class LostNum(NodeNum):
    """
        Type for lacking leaves.
    """

    @property
    def num(self):
        return self._num

    @num.setter
    def num(self, value):
        self._num = None

    @staticmethod
    def _get_num(other):
        return None
    
    @property
    def name(self):
        return "l"

    @property
    def is_last(self):
        return True

    def __bool__(self):
        return False


class NodeContacts:
    """
    Class for node's representation in nodes dictionary.
    """
    def __init__(self,
                 parent: QubitNum | str = 0,
                 childs: list[QubitNum | BranchNum | LostNum] = None
                 ):
        if parent == ROOT:
            self.parent = parent
        else:
            self.parent = QubitNum(parent)
        if childs is None:
            self._childs = []
        else:
            self._childs = childs
        self.childs = self._childs

    @property
    def childs(self):
        return self._childs

    def child_index(self, child):
        return self.childs.index(child)

    @staticmethod
    def obj_to_child(value):
        if isinstance(value, int):
            return QubitNum(value)
        elif isinstance(value, NodeNum):
            return value
        raise Exception("Unknown object for childs. Should be int or NodeNum")

    @childs.setter
    def childs(self, value):
        childs = [None, None, None]
        i = 0
        for child in value:
            childs[i] = self.obj_to_child(child)
            i += 1
        for j in range(i, 3):
            childs[j] = LostNum()
        self._childs = tuple(childs)

    def __setitem__(self, key, val):
        childs = list(self.childs)
        childs[key] = self.obj_to_child(val)
        self.childs = tuple(childs)
        return val

    def __getitem__(self, key):
        return self.childs[key]

    def __str__(self):
        return str((self.parent,) + self.childs)

    def __repr__(self):
        return str(self.parent) + "," + str(self.childs)
