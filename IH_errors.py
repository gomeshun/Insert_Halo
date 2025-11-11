
class StatusException(BaseException):
    """
    Is raised when the global status of Insert Halo indicates that some process is still running while another process is being called.
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

