class SkyCatalogsException(Exception):
    pass

class NoSchemaVersionError(SkyCatalogsException):
    def __init__(self, msg):

        if not msg:
            msg = 'Schema version unspecified'
        self.msg = msg
        super().__init__(self.msg)
