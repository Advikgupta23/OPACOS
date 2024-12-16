import sys

galactic_long = sys.argv[1]
galactic_lat = sys.argv[2]
survey_area = sys.argv[3]
with open('/Users/advik/GalaxiaData/Examples/myparameterfile','r') as file:
	lines = file.readlines()
	lines[11] = f"longitude                           {galactic_long}\n"
	lines[12] = f"latitude                            {galactic_lat}\n"
	lines[13] = f"surveyArea                          {survey_area}\n"
with open('/Users/advik/GalaxiaData/Examples/myparameterfile', 'w') as file:
	file.writelines(lines)
