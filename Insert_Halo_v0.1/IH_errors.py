
class StatusException(BaseException):
    """
    Is raised when the global status of Insert Halo indicates that some process is still running while another process is being called.
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

class StageException(BaseException):
    """
    Is raised when the current function of Insert Halo requires some previous computation to have been run first, which has not happened.
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

