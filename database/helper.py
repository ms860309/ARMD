import datetime
import sys
import platform


def highlight_text(text, width=80, level=0):
    """
    This code is from pysisyphus.
    https://github.com/eljost/pysisyphus
    """
    levels = {
        #  horizontal
        #        vertical
        0: ("#", "#"),
        1: ("-", "|"),
    }
    full_length = len(text) + 4
    pad_len = width - full_length
    pad_len = (pad_len - (pad_len % 2)) // 2
    pad = " " * pad_len
    hchar, vchar = levels[level]
    full_row = hchar * full_length
    highlight = (
        f"""{pad}{full_row}\n{pad}{vchar} {text.upper()} {vchar}\n{pad}{full_row}"""
    )
    return highlight


def print_header():
    """
    Generated from Text to ASCII Art Generator (TAAG)
    """
    logo = """
    ___         __                        __           __                                  
   /   | __  __/ /_____  ____ ___  ____ _/ /____  ____/ /                                  
  / /| |/ / / / __/ __ \/ __ `__ \/ __ `/ __/ _ \/ __  /                                   
 / ___ / /_/ / /_/ /_/ / / / / / / /_/ / /_/  __/ /_/ /                                    
/_/  |_\__,_/\__/\____/_/_/_/_/_/\__,_/\__/\___/\__,_/      _                              
                        /  |/  /__  _____/ /_  ____ _____  (_)________ ___                 
                       / /|_/ / _ \/ ___/ __ \/ __ `/ __ \/ / ___/ __ `__ \                
                      / /  / /  __/ /__/ / / / /_/ / / / / (__  ) / / / / /                
                     /_/  /_/\___/\___/_/ /_/\__,_/_/ /_/_/____/_/ /_/ /_/                 
                                               ____  _                                     
                                              / __ \(_)_____________ _   _____  _______  __
                                             / / / / / ___/ ___/ __ \ | / / _ \/ ___/ / / /
                                            / /_/ / (__  ) /__/ /_/ / |/ /  __/ /  / /_/ / 
                                           /_____/_/____/\___/\____/|___/\___/_/   \__, /  
                                                                                  /____/   """

    logo_2 = """
   ___       __                  __         __               
  / _ |__ __/ /____  __ _  ___ _/ /____ ___/ /               
 / __ / // / __/ _ \/  ' \/ _ `/ __/ -_) _  /                
/_/ |_\_,_/\__/\___/_/_/_/\_,_/\__/\__/\_,_/_                
             /  |/  /__ ____/ /  ___ ____  (_)__ __ _        
            / /|_/ / -_) __/ _ \/ _ `/ _ \/ (_-</  ' \       
           /_/  /_/\__/\__/_//_/\_,_/_//_/_/___/_/_/_/       
                        / _ \(_)__ _______ _  _____ ______ __
                       / // / (_-</ __/ _ \ |/ / -_) __/ // /
                      /____/_/___/\__/\___/___/\__/_/  \_, / 
                                                      /___/  
                                                            """
    vi = sys.version_info
    sv = f"{vi.major}.{vi.minor}.{vi.micro}"  # Python

    print(
        f"{logo_2}\n\nPython {sv}\n"
        f"Executed at {datetime.datetime.now().strftime('%c')} on '{platform.node()}'\n"
        f"Platform: {platform.platform()}\n"
    )
