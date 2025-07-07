class CuPyManager:
    _instance = None

    def __init__(self):
        self._cp = None
        self._HAS_CUPY = False

    def _load_cupy(self):
        if self._cp is None:  # Only attempt load once
            try:
                import cupy as cp

                self._cp = cp
                self._HAS_CUPY = True
            except (ImportError, OSError):
                self._cp = None
                self._HAS_CUPY = False

    @property
    def cp(self):
        self._load_cupy()
        return self._cp

    @property
    def HAS_CUPY(self):
        self._load_cupy()
        return self._HAS_CUPY


# Singleton instance
cupy_manager = CuPyManager()
