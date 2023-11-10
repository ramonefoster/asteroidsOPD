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

def is_numeric(input):
    """
    Check if the input is a numeric value.
    Parameters:
        input (int or float): The value to be checked.
    Returns:
        bool: True if the input is numeric (int or float), False otherwise.
    """
    if isinstance(input, (int, float)):
        return True
    else:
        return False
    
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

def hours_to_hms(hours, decimal_digits=0):
    """
    Converts Float Hour to string Hour, in format hh:mm:ss:cc
    :param hours: Hours (float)
    """
    if is_numeric(hours):        
        sign = "-" if hours < 0 else ""
        hours = abs(hours)
        whole_hours = int(hours)
        fractional_hours = hours - whole_hours

        minutes = int(fractional_hours * 60)
        fractional_minutes = fractional_hours * 60 - minutes

        seconds = int(fractional_minutes * 60)
        fractional_seconds = fractional_minutes * 60 - seconds

        seconds_str = f"{seconds:02}.{int(fractional_seconds * (10 ** decimal_digits)):02d}"

        time_string = f"{sign}{whole_hours:02}:{minutes:02}:{seconds_str}"
        
        return time_string
    else:
        return None

def degrees_to_dms(degrees):
    """
    Converts Degrees to string, in format dd:mm:ss:cc
    :param hours: Degrees (float)
    """
    if is_numeric(degrees):
        sign = "-" if degrees < 0 else "+"
        degrees = abs(degrees)
        degrees_int = int(degrees)
        minutes = int((degrees - degrees_int) * 60)
        seconds = int(((degrees - degrees_int) * 60 - minutes) * 60)
        seconds_decimal = int((((degrees - degrees_int) * 60 - minutes) * 60 - seconds) * 100)

        # Formated value
        degrees_string = f'{sign}{degrees_int:02}:{minutes:02}:{seconds:02}.{seconds_decimal:02}'

        return degrees_string
    else:
        return None