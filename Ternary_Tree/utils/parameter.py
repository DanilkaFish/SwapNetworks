class Parameter:
    """
        var = coef * name
    """
    def __init__(self, name: str='t'):
        self.name = name
        self.coef = 1
        