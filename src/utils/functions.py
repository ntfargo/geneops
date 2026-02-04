import os

def lower_list(input: list): return [v.lower() for v in input]

def lower_dict(input: dict): 
    return dict((k.lower(), v.lower()) for k, v in input.items())

def find_test_dir(start_dir: str | None = None) -> str:
    if not start_dir:
        start_dir = "."

    target = os.path.abspath(start_dir)
    while True:
        if os.path.isdir(os.path.join(target, "src")) and os.path.isdir(
            os.path.join(target, "tests")
        ):
            return os.path.abspath(os.path.join(target, "tests"))
        new, tmp = os.path.split(target)
        if target == new:
            break
        target = new
    raise ValueError(
        f"Could not find \"tests\" directory starting from {start_dir}"
    )