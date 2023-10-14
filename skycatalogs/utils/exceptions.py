__all__ = ['SkyCatalogsException', 'NoSchemaVersionError',
           'ConfigDuplicateKeyError']


class SkyCatalogsException(Exception):
    pass


class NoSchemaVersionError(SkyCatalogsException):
    def __init__(self, msg):

        if not msg:
            msg = 'Schema version unspecified'
        self.msg = msg
        super().__init__(self.msg)


class ConfigDuplicateKeyError(SkyCatalogsException):
    def __init__(self, key):
        self.msg = f'Cannot add duplicate key {key} to config'
        super().__init__(self.msg)


class SkyCatalogsRuntimeError(SkyCatalogsException):
    def __init__(self, msg):
        if not msg:
            msg = 'skyCatalogs runtime error'
        self.msg = msg
        super().__init__(self.msg)
