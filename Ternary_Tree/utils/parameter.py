class Parameter:
    """
        var = coef * name
    """
    def __init__(self, name: str='t'):
        self.name = name
        self.coef = 1

    def __repr__(self):
        return self.name
        
    def __str__(self):
        return self.name