def check_format(input):
    """
    Check the format of a given input and extract components.

    The expected format is 'x1:x2:x3' or 'x1 x2 x3',
    where 'x1', 'x2', and 'x3' are separated by a colon (':') or a space (' ').

    Parameters:
        input (str): The input to be checked and parsed.

    Returns:
        False: If the input does not match the expected format.
        list: A list containing the three components extracted from the input if the format is correct.
    """
    separators = [':', ' ']
    separator = None
    for sep in separators:
        if sep in input:
            separator = sep
            break

    if separator is None:
        return False

    components = input.split(separator)

    # Check for correct format
    if len(components) != 3:
        return False
    else:
        return components

def hms_to_hours(time_string):
    """
    Converts Hours string to float
    :param time_string: Hours String (hh:mm:ss.ss)
    """        
    # Verify separator
    components = check_format(time_string)

    if components:
        hours = abs(int(components[0]))
        minutes = int(components[1])
        seconds = float(components[2])

        total_hours = hours + minutes / 60 + seconds / 3600

        sign = -1 if "-" in time_string else 1
        return sign*total_hours
    else:
        return None

def dms_to_degrees(degrees_string):
    """
    Converts Degrees string to float
    :param degrees_string: Degrees String (dd:mm:ss.ss)
    """
    # Verify separator
    components = check_format(degrees_string)

    if components:
        degrees_int = abs(int(components[0]))
        minutes = int(components[1])    
        seconds = float(components[2])

        degrees = degrees_int + minutes / 60 + seconds / 3600

        sign = -1 if "-" in degrees_string else 1
        return sign*degrees
    else:
        return None