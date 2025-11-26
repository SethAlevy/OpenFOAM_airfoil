from datetime import datetime


class SimpleLogger:
    @staticmethod
    def log(message: str) -> None:
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {message}")

    @staticmethod
    def warning(message: str) -> None:
        # ANSI escape code for yellow
        YELLOW = "\033[93m"
        RESET = "\033[0m"
        print(
            f"{YELLOW}[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] WARNING: \
                  {message}{RESET}")
