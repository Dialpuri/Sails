import os
import tempfile

def swap_files():
    file1 = "_pyproject.toml"
    file2 = "pyproject.toml"

    file3 = "setup.py"
    file4 = "_setup.py"


# Check if both files exist
    if not os.path.isfile(file1) or not os.path.isfile(file2):
        return

    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
        temp_name = tmpfile.name

    try:
        os.rename(file1, temp_name)
        os.rename(file2, file1)
        os.rename(temp_name, file2)

        if os.path.exists(file3):
            os.rename(file3, file4)
        else:
            os.rename(file4, file4)

    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        # Ensure that the temporary file is removed in case of failure
        if os.path.exists(temp_name):
            os.remove(temp_name)

if __name__ == "__main__":
    swap_files()
