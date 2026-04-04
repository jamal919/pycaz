# -*- coding: utf-8 -*-
import logging
import signal
import time
from functools import wraps


def retry(func: callable, retries: int = 5, delay: int = 1, exceptions: Exception = (Exception,),
          logger: logging.Logger = None):
    """Helper function to retry a function

    Args:
        func (callable): A callable function
        retries (int, optional): Number of retries. Defaults to 3.
        delay (int, optional): Delay before next try. Defaults to 1.
        exceptions (Exception, optional): Which exceptions to consider. Defaults to (Exception,).
        logger (logging.Logger, optional): Data logger. Defaults to None.
    """
    if logger is None:
        logger = logging.getLogger("retry")
    else:
        logger = logger

    for attempt in range(1, retries + 1):
        try:
            return func()
        except exceptions as e:
            logger.info(f"Attempt {attempt}/{retries} failed with {e}")
            if attempt == retries:
                raise Exception("All retries failed")
            time.sleep(delay)
            continue


def timeout(seconds: int):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Define the handler for the alarm signal
            def handler(signum, frame):
                raise TimeoutError(f"Function '{func.__name__}' timed out after {seconds} seconds.")

            # Set the signal handler and the alarm
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(seconds)

            try:
                result = func(*args, **kwargs)
            finally:
                # Disable the alarm regardless of success or failure
                signal.alarm(0)
            return result

        return wrapper

    return decorator
