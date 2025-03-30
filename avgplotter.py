"""
	====================================================================================
	Plots the frequency density of the dihedral (torsion/Ramachandran) angles from the 
	Top8000 peptide database as a smoothed background histogram, along with two contour 
	lines.

	The averaged dihedral angles (phi, psi) from multiple user-supplied protein structures 
	are then overlaid as a scatter plot. In addition, the average angles for each PDB 
	are saved to a text file (and optionally as a CSV).

	Colours and recommended parameters can be easily adjusted towards the end of the 
	script.
	
	Version 2.1.0:
	 - Now processes multiple PDB files (comma separated) and plots the averaged angles.
	 - Retains previous functionality with a few modifications.
	
	Author information:
	 - Joseph I. J. Ellaway
	 - josephellaway@gmail.com
	 - https://github.com/Joseph-Ellaway
	====================================================================================
"""

import os
import matplotlib.pyplot as plt
import pandas as pd
# Package functions
from DihedralCalculator import *
from PlotterFunctions import *
from RamaArgumentParser import *

# Main function
def main(pdb, itmod, model_num, itchain, chain_num, plot_type, out_dir, verb, save, file_type):
    ########################################################
    #           PRE-PROCESSING USER INPUT                #
    ########################################################
    
    # Expecting pdb to be a comma-separated list if multiple files are provided.
    pdb_list = [p.strip() for p in pdb]

    # Available plot options
    options = ["All", "General", "Glycine", "Proline", "Pre-proline", "Ile-Val"]
    plot_type = options[int(plot_type)]
    
    # Set a base name for the plot (note: not tied to a single pdb id)
    plot_name = os.path.join(out_dir, 'averaged_' + plot_type + "RamachandranPlot_tmp")
    
    # Prepare lists to store averaged phi/psi angles and pdb ids
    averaged_phi = []
    averaged_psi = []
    pdb_ids = []
    
    ########################################################
    #           PROCESS EACH PDB FILE                    #
    ########################################################
    
    for pdb_file in pdb_list:
        VerboseStatement(verb, "Importing " + pdb_file)
        
        # Extract dihedral angles for this pdb file
        userpdb_df = ExtractDihedrals(
            pdb_file_name=pdb_file, 
            iter_models=itmod, 
            model_number=model_num, 
            iter_chains=itchain, 
            chain_id=chain_num
        )
        userpdb_df = userpdb_df.dropna()
        
        # Filter by residue type if needed
        if plot_type != "All":
            userpdb_df = userpdb_df.loc[userpdb_df["type"] == plot_type]
        
        if userpdb_df.empty:
            VerboseStatement(verb, f"No valid dihedral data found for {pdb_file}. Skipping.")
            continue
        
        # Calculate the arithmetic averages (note: for circular data you may need a different approach)
        avg_phi = userpdb_df["phi"].mean()
        avg_psi = userpdb_df["psi"].mean()
        
        averaged_phi.append(avg_phi)
        averaged_psi.append(avg_psi)
        
        # Store the pdb id (using the file name without extension)
        pdb_id = os.path.basename(pdb_file).split('.')[0]
        pdb_ids.append(pdb_id)
    
    if not averaged_phi:
        print("No valid dihedral data found for any pdb file. Exiting.")
        return
    
    VerboseStatement(verb, "Dihedral angles calculated for all pdb files")
    
    # Optionally save the averaged values as a CSV file
    if save:
        averages_csv = os.path.join(out_dir, 'averaged_angles.csv')
        df_avg = pd.DataFrame({
            "pdb_id": pdb_ids,
            "avg_phi": averaged_phi,
            "avg_psi": averaged_psi
        })
        VerboseStatement(verb, "Saving averaged angles CSV as: " + averages_csv)
        df_avg.to_csv(averages_csv, index=False)
    
    ########################################################
    #           IMPORTING REFERENCE DATA                 #
    ########################################################
    
    VerboseStatement(verb, "Importing Top8000 library")
    top8000_df = SelectAngles(
        pd.read_csv("Top8000_DihedralAngles.csv.gz", compression="gzip"), 
        plot_type
    )
    
    ########################################################
    #           GENERATE BACKGROUND IMAGE              #
    ########################################################
    
    # Recommended adjustable parameters for the background
    background_colour = "Blues"  # See: https://matplotlib.org/stable/tutorials/colors/colormaps.html
    
    VerboseStatement(verb, "Generating background of favoured regions")
    MakeBackground(top8000_df, plot_type, plot_name, background_colour)
    
    ########################################################
    #           PLOTTING THE DATA                        #
    ########################################################
    
    # Plot parameters
    figure_size = (5, 5)             # Figure size in inches
    contour_level_inner = 96         # Percentile for inner contour
    contour_level_outer = 15         # Percentile for outer contour
    contour_line_color_inner = "#DFF8FB"
    contour_line_color_outer = "#045E93"
    out_resolution = 96              # Output figure resolution
    data_point_colour = "#D4AB2D"      # Colour for averaged points
    data_point_edge_colour = "#3c3c3c"
    
    VerboseStatement(verb, "Plotting averaged Ramachandran diagram")
    plt.style.use("seaborn-v0_8-poster")
    fig, ax = plt.subplots(1, 1, figsize=figure_size, tight_layout=True)
    
    # ADDING CONTROURS
    AddContour(ax, top8000_df, contour_level=contour_level_inner, line_colour=contour_line_color_inner)
    AddContour(ax, top8000_df, contour_level=contour_level_outer, line_colour=contour_line_color_outer, contour_alpha=0.3)
    
    # ADDING BACKGROUND IMAGE (generated by MakeBackground)
    ax.imshow(plt.imread(str(plot_name + ".png")), extent=[-195, 195, -195, 195], zorder=1)
    
    # ADD GRIDLINES & FORMAT AXES
    AddGridLines(ax)
    FormatAxis(ax)
    
    # PLOTTING THE AVERAGED DIHEDRAL ANGLES FOR EACH PDB
    ax.scatter(averaged_phi, averaged_psi, s=15, color=data_point_colour, 
               zorder=4, linewidths=0.5, edgecolor=data_point_edge_colour)
    
    # SAVE THE PLOT
    VerboseStatement(verb, "Saving plot")
    if file_type == "png":
        final_plot_file = str(plot_name[:-4] + ".png")
        plt.savefig(final_plot_file, dpi=out_resolution, bbox_inches=0, pad_inches=None)
    else:
        final_plot_file = str(plot_name[:-4] + '.' + file_type)
        plt.savefig(final_plot_file, bbox_inches=0, pad_inches=None)
    
    # Remove temporary background image
    os.remove(str(plot_name + '.png'))
    plt.close()
    
    ########################################################
    #           OUTPUT THE AVERAGE ANGLES              #
    ########################################################
    
    averages_file = str(plot_name[:-4] + '_averages.txt')
    with open(averages_file, "w") as file:
        for pdb_id, avg_phi_val, avg_psi_val in zip(pdb_ids, averaged_phi, averaged_psi):
            file.write("PDB: {} | Average phi: {} | Average psi: {}\n".format(pdb_id, avg_phi_val, avg_psi_val))
    
    print("Averages saved to", averages_file)
    print("Done. \nRamachandran plot saved to", final_plot_file)


########################################################
#           RUNNING SCRIPT                           #
########################################################

if __name__ == "__main__":
    if __name__ == "__main__":
        args = CollctUserArgs()
    main(
        args.pdb,          # list of PDB files
        args.models,       # number of models to iterate through
        1,                 # model_num (default, change as needed)
        1,                 # itchain (default, change as needed)
        args.chains,       # chain identifier
        args.plot_type,    # plot type (integer corresponding to type)
        args.out_dir,      # output directory
        args.verbose,      # verbose flag
        args.save,         # save flag
        args.file_type     # output file type (png/pdf)
    )

else:
    pass


