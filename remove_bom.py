# coding: utf-8
# -----------------
# created: 30.05.2024
# author: Robby Heusequin
# purpose: Remove BOM encoding from a GPX file
# -----------------
# changes: 
#
# ----------------------------------------------

def remove_bom_encoding(file_path):
    with open(file_path, 'rb') as f:
        gpx_content = f.read()
    # If BOM is present, remove it
    if gpx_content.startswith(b'\xef\xbb\xbf'):
        gpx_content = gpx_content[3:]
    with open(file_path, 'wb') as f:
        f.write(gpx_content)

file_path = "ENTER GPX FILE"
remove_bom_encoding(file_path)