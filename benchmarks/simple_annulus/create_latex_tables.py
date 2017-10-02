import os
import re
import numpy
grid_resolutions=['2', '3', '4', '5', '6', '7']
interpolation_schemes = ['bilinear', 'cell_average']
density_directory_names=['rho_0', 'rho_1', 'rho_r', 'rho_r^2']
densities=['rho = 0', 'rho = 1', 'rho = r', 'rho = r\\^2']
FEM_types=['Q2_Q1', 'Q2_P-1', 'Q3_Q2', 'Q3_P-2', 'Q1_P0']
FEM_names=['Q2xQ1', 'Q2xP-1', 'Q3xQ2', 'Q3xP-2', 'Q1xP0']

#Some sources of information used for creating this script:
#https://stackoverflow.com/questions/19859282/check-if-a-string-contains-a-number
#https://stackoverflow.com/questions/23636509/python-convert-string-in-scientific-notation-to-float
#https://stackoverflow.com/questions/6913532/display-a-decimal-in-scientific-notation
#https://stackoverflow.com/questions/15263597/convert-floating-point-number-to-certain-precision-then-copy-to-string
def main():
		
	latex_tables_document = []
	document_start = readTemplateFile()
	latex_tables_document.append(document_start)
	
	latex_tables = construct_direct_method_tables()
	latex_tables_document.append(latex_tables)
	particle_error_tables = construct_particle_error_tables()
	latex_tables_document.append(particle_error_tables)		
	write_tables_to_tex_file(latex_tables_document)
		
def write_tables_to_tex_file(latex_tables_document):
	file_writer = open("error_tables.tex", "w")
	
	for line in latex_tables_document[0]:
		file_writer.write(line)
		
	#Write direct error tables to the latex file
	for table in latex_tables_document[1]:
		for line in table:
			file_writer.write(line)
			file_writer.write("\n")
			
	#Write particle error tables to the latex file
	for table in latex_tables_document[2]:
		for line in table:
			file_writer.write(line)
			file_writer.write("\n")
		
	file_writer.write("\\end{document}")
	file_writer.write("\n")
	file_writer.close()
	

	
	
	
def construct_particle_error_tables():
	latex_tables_document = []

	for i in range(0, len(density_directory_names)):
		density_dir = density_directory_names[i]
		#Collect the errors
		FEM_interpolation_pressure_errors = []
		FEM_interpolation_velocity_errors = []
				
		for FEM_type in FEM_types:
			interpolation_pressure_errors = []
			interpolation_velocity_errors = []
			for interpolation_scheme in interpolation_schemes:
				
				pressure_errors = []
				velocity_errors = []
				for resolution in grid_resolutions:
					parameter_file_name = "simple_annulus_" + resolution + "x" + resolution + ".out"
					parameter_file_path = os.path.join(os.getcwd(), "compositional_fields", density_dir, FEM_type, interpolation_scheme, resolution, parameter_file_name)
					
					error_data = readOutputFile(parameter_file_path)
				
					#Get the L2 error data for the last time step if there is any error data
					if len(error_data) == 0:
						velocity_errors.append("")
						pressure_errors.append("")
					else:
						
						#print(error_data[-2:0])
						velocity_errors.append(error_data[-1][-2])
						pressure_errors.append(error_data[-1][-1])
						
				interpolation_pressure_errors.append(pressure_errors)	
				interpolation_velocity_errors.append(velocity_errors)	
			


			FEM_interpolation_pressure_errors.append(interpolation_pressure_errors)
			FEM_interpolation_velocity_errors.append(interpolation_velocity_errors)
			
			
		#Generate tables from the errors
		for FEM in range(0, len(FEM_types)):	
			velocity_error_table = construct_particle_method_table_latex("velocity_errors", FEM_interpolation_velocity_errors[FEM], densities[i], FEM_names[FEM])
			pressure_error_table = construct_particle_method_table_latex("pressure_errors", FEM_interpolation_pressure_errors[FEM], densities[i], FEM_names[FEM])
			latex_tables_document.append(velocity_error_table)
			latex_tables_document.append(pressure_error_table)
		
	return latex_tables_document
	
def construct_direct_method_tables():
	latex_tables_document = []

	for i in range(0, len(density_directory_names)):
		density_dir = density_directory_names[i]
		#Collect the errors
		FEM_pressure_errors = []
		FEM_velocity_errors = []
			
		for FEM_type in FEM_types:
			pressure_errors = []
			velocity_errors = []
			for resolution in grid_resolutions:
				parameter_file_name = "simple_annulus_" + resolution + "x" + resolution + ".out"
				parameter_file_path = os.path.join(os.getcwd(), density_dir, FEM_type, resolution, parameter_file_name)
				error_data = readOutputFile(parameter_file_path)

				#Get the L2 error data for the last time step if there is any error data
				if len(error_data) == 0:
					velocity_errors.append("")
					pressure_errors.append("")
				else:
					velocity_errors.append(error_data[-1][-2])
					pressure_errors.append(error_data[-1][-1])
					
			FEM_pressure_errors.append(pressure_errors)	
			FEM_velocity_errors.append(velocity_errors)	
			
		velocity_error_table = construct_direct_method_table_latex("velocity_errors", FEM_velocity_errors, densities[i])
		pressure_error_table = construct_direct_method_table_latex("pressure_errors", FEM_pressure_errors, densities[i])
		latex_tables_document.append(velocity_error_table)
		latex_tables_document.append(pressure_error_table)
	
	
	return latex_tables_document
	
	
	
	
#Build a latex table to display ASPECT errors collected using the direct method
#table_type is either "velocity_errors" or "pressure_errors"
#aspect_errors is a 2d array containing the error data (pressure or velocity) for all finite
# element types and grid resolutions 

#Returns an array, where each entry in the array is a line of latex code
#for the latex table
def construct_direct_method_table_latex(table_type, aspect_errors, density):
	latex_table = []
	#Generate the beginning of the latex table
	latex_table.append("Density: " + density + "\\newline")
	if table_type == "velocity_errors":
		latex_table.append("%%%%%%%%%   u_L2")
	else:
		latex_table.append("%%%%%%%%%   p_L2")
		
	latex_table.append("\\begin{table}[ht!]")
	latex_table.append("\\centering")
	
	if table_type == "velocity_errors":
		latex_table.append("\\caption{Velocity Errors}")
	else:
		latex_table.append("\\caption{Pressure Errors}")
		
	latex_table.append("\\vskip 6pt")
	
	if table_type == "velocity_errors":
		latex_table.append("\\label{Velocity Errors}")
	else:
		latex_table.append("\\label{Pressure Errors}")
	
	latex_table.append("\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}")
	latex_table.append("\\hline")
	
	if table_type == "velocity_errors":
		latex_table.append("\\multicolumn{9}{|c|}{Velocity Errors for Direct Method}  \\\ \hline ")
	else:
		latex_table.append("\\multicolumn{9}{|c|}{Pressure Errors for Direct Method} \\\ \hline")
	
	latex_table.append("   h & \\multicolumn{2}{c|}{$\Q_{2}\mathrm{x}\Q_{1}$} & \multicolumn{2}{c|}{$\Q_{2}\mathrm{x}P_{-1}$}& \multicolumn{2}{c|}{$\Q_{3}\mathrm{x}\Q_{2}$}& \multicolumn{2}{c|}{$\Q_{3}\mathrm{x}P_{-2}$} & \multicolumn{2}{c|}{$\Q_{1}\mathrm{x}P_{0}$} \\\ \hline")
	latex_table.append("$h$ &  error          & rate           & error          & rate           & error          & rate      & error  & rate  & error  & rate \\\ \hline")
	#Add the error data to the table
	for res in range(0, len(grid_resolutions)):
		latex_table_row = grid_resolutions[res] + "   & "
		for fem in range(0, len(FEM_types)):	
			latex_table_row = str(latex_table_row) + aspect_errors[fem][res] + " & "
			#Add in the convergence rate if applicable
			if res == 0 or aspect_errors[fem][res - 1] == "" or aspect_errors[fem][res] == "":
				latex_table_row = latex_table_row + "x.xx "
			else:
				rate = numpy.log2(float(aspect_errors[fem][res - 1])/float(aspect_errors[fem][res]))
				latex_table_row = latex_table_row + str(round(rate,1)) + "0"
				
			if fem == len(FEM_types) - 1:
				latex_table_row = latex_table_row + " \\\ \hline"
			else:
				latex_table_row = latex_table_row + " & "
		latex_table.append(latex_table_row)
		
	latex_table.append("\end{tabular}")
	latex_table.append("\end{table}")
	latex_table.append("")
	
	return latex_table
	
def construct_particle_method_table_latex(table_type, aspect_errors, density, FEM):
	latex_table = []
	#Generate the beginning of the latex table
	latex_table.append("Density: " + density + "\\newline")
	if table_type == "velocity_errors":
		latex_table.append("%%%%%%%%%   u_L2")
	else:
		latex_table.append("%%%%%%%%%   p_L2")
		
	latex_table.append("\\begin{table}[ht!]")
	latex_table.append("\\centering")
	
	if table_type == "velocity_errors":
		latex_table.append("\\caption{Velocity Errors}")
	else:
		latex_table.append("\\caption{Pressure Errors}")
		
	latex_table.append("\\vskip 6pt")
	
	if table_type == "velocity_errors":
		latex_table.append("\\label{Table:" + FEM + ": Velocity Errors}")
	else:
		latex_table.append("\\label{Table:" + FEM +  " Pressure Errors}")
	
	latex_table.append("\\begin{tabular}{|c|c|c|c|c|}")
	latex_table.append("\hline")
	
	if table_type == "velocity_errors":
		latex_table.append("\\multicolumn{5}{|c|}{" + FEM + ": Velocity Errors for Particle Method}  \\\ \hline ")
	else:
		latex_table.append("\\multicolumn{5}{|c|}{" + FEM + ": Pressure Errors for Particle Method}  \\\ \hline ")
	
	latex_table.append("    $h$  & \multicolumn{2}{c|}{cell average} & \multicolumn{2}{c|}{bilinear} \\\ \hline")
	latex_table.append("$h$ &  error          & rate           & error          & rate  \\\ \hline")
	#Add the error data to the table
	for res in range(0, len(grid_resolutions)):
		latex_table_row = grid_resolutions[res] + "   & "
		for interpolation in range(0, len(interpolation_schemes)):	
			latex_table_row = str(latex_table_row) + aspect_errors[interpolation][res] + " & "
			#Add in the convergence rate if applicable
			if res == 0 or aspect_errors[interpolation][res - 1] == "" or aspect_errors[interpolation][res] == "":
				latex_table_row = latex_table_row + "x.xx "
			else:
				rate = numpy.log2(float(aspect_errors[interpolation][res - 1])/float(aspect_errors[interpolation][res]))
				latex_table_row = latex_table_row + str(round(rate,1)) + "0"
				
			if interpolation == len(interpolation_schemes) - 1:
				latex_table_row = latex_table_row + " \\\ \hline"
			else:
				latex_table_row = latex_table_row + " & "
		latex_table.append(latex_table_row)
		
	latex_table.append("\end{tabular}")
	latex_table.append("\end{table}")
	latex_table.append("")
	
	return latex_table

#Read the given output file, and return an array of arrays containing all of the error
#data in the file.
#Each entry in the first array contains all of the error data for a single time step
#Each sub array contains the individual errors (u_L1, p_L1, u_L2, p_L2)
def readOutputFile(file_name):
	file_info = []
	if os.path.isfile(file_name):
		file_reader = open(file_name, "r")
		
		line = ""
		
		for line in file_reader:
			if "Errors" in line:
				# Split the error output on two separators, "," and " ".
				parsed_error_data = re.split("[,\s]+", line.strip())
				#The last four items in the split string are the pressure and velocity errors
				errors = parsed_error_data[-4:]
				file_info.append(errors)
	
	return file_info
	
# Read the EGP Latex template file
def readTemplateFile():
	file_reader = open("Particle_Interpolation_Time_Dependent_Flow.tex", "r")
	file_info = []
	
	for line in file_reader:
		file_info.append(line)
		
	return file_info






























if __name__ == '__main__':
	main()