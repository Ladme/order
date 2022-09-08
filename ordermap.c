// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

// Creates lipid order map for every bond of every requested residue.

#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "general.h"

const char VERSION[] = "v2022/08/26";

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;
// inverse size of a grid tile for order map calculation
static const int GRID_TILE = 10;



/* Get user-defined grid dimension. */
int get_range(const char *optarg, float array[2])
{
    if (sscanf(optarg, "%f-%f", &array[0], &array[1]) != 2 &&
        sscanf(optarg, "%f - %f", &array[0], &array[1]) != 2 &&
        sscanf(optarg, "%f %f", &array[0], &array[1]) != 2) {
        fprintf(stderr, "Could not understand grid dimension specifier.\n");
        return 1;
    }

    return 0;
}


/*
 * Parses command line arguments.
 * Returns zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments(
        int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **itp_file,
        char **atoms,
        char **phosphates,
        float *array_dimx,
        float *array_dimy,
        char **output_directory,
        int *nan_limit) 
{
    int gro_specified = 0, xtc_specified = 0;

    int opt = 0;
    while((opt = getopt(argc, argv, "c:f:n:i:a:p:x:y:o:l:h")) != -1) {
        switch (opt) {
        // help
        case 'h':
            return 1;
        // gro file to read
        case 'c':
            *gro_file = optarg;
            gro_specified = 1;
            break;
        // xtc file to read
        case 'f':
            *xtc_file = optarg;
            xtc_specified = 1;
            break;
        // ndx file to read
        case 'n':
            *ndx_file = optarg;
            break;
        // itp file to use
        case 'i':
            *itp_file = optarg;
            break;
        // atoms to analyze
        case 'a':
            *atoms = optarg;
            break;
        // phosphates identifier
        case 'p':
            *phosphates = optarg;
            break;
        // grid dimensions
        case 'x':
            if (get_range(optarg, array_dimx) != 0) return 1;
            break;
        case 'y':
            if (get_range(optarg, array_dimy) != 0) return 1;
            break;
        // output directory
        case 'o':
            *output_directory = optarg;
            break;
        // nan limit
        case 'l':
            sscanf(optarg, "%d", nan_limit);
            break;
        default:
            return 1;
        }
    }

    if (!gro_specified || !xtc_specified) {
        fprintf(stderr, "Gro file and xtc file must always be supplied.\n");
        return 1;
    }
    return 0;
}

void print_usage(const char *program_name)
{
    printf("Usage: %s -c GRO_FILE -f XTC_FILE [OPTION]...\n", program_name);
    printf("\nOPTIONS\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-i STRING        itp file to read (default: martini_v3.0.0_phospholipids_v1.itp)\n");
    printf("-o STRING        output directory (default: order_maps)\n");
    printf("-a STRING        selection of atoms to be used for analysis (default: Membrane)\n");
    printf("-p STRING        identifier of lipid heads (default: name PO4)\n");
    printf("-x FLOAT         grid dimension in x axis (default: box size from gro file)\n");
    printf("-y FLOAT         grid dimension in y axis (default: box size from gro file)\n");
    printf("-l INTEGER       NAN limit: how many bonds must be detected in a grid tile\n");
    printf("                 to calculate lipid order for this tile (default: 30)\n");
    printf("\n");
}

/*
 * Prints parameters that the program will use for the calculation.
 */
void print_arguments(
        FILE *stream,
        const char *gro_file,
        const char *xtc_file,
        const char *ndx_file,
        const char *itp_file,
        const char *atoms,
        const char *phosphates,
        const float array_dimx[2],
        const float array_dimy[2],
        const char *output_directory,
        const int nan_limit) 
{
    fprintf(stream, "\nParameters for Lipid Order Map calculation:\n");
    fprintf(stream, ">>> gro file:         %s\n", gro_file);
    if (xtc_file == NULL) fprintf(stream, ">>> xtc file:         ----\n");
    else fprintf(stream, ">>> xtc file:         %s\n", xtc_file);
    fprintf(stream, ">>> ndx file:         %s\n", ndx_file);
    fprintf(stream, ">>> itp file:         %s\n", itp_file);
    fprintf(stream, ">>> output directory: %s\n", output_directory);
    fprintf(stream, ">>> atoms:            %s\n", atoms);
    fprintf(stream, ">>> phosphates:       %s\n", phosphates);
    fprintf(stream, ">>> grid dimensions:  x: %.1f - %.1f, y: %.1f - %.1f\n", 
            array_dimx[0], array_dimx[1], array_dimy[0], array_dimy[1]);
    fprintf(stream, ">>> NAN limit:        %d\n", nan_limit);
    fprintf(stream, "\n");
}

/*
 * Converts index of an array to coordinate.
 */
static inline float index2coor(int x, float minx)
{
    return (float) x / GRID_TILE + minx;
}
        
/*
 * Converts coordinate to an index array.
 */
static inline size_t coor2index(float x, float minx)
{
    return (size_t) roundf((x - minx) * GRID_TILE);
}


/* Creates an order or samples map for all residue types and bonds. 
 * maptype: 0 -> order map, 1 -> samples map
 */
dict_t *create_map(
        const list_t *resnames, 
        const dict_t *bonds,
        const size_t n_tiles,
        const int maptype)
{
    dict_t *map_dict = dict_create();

    for (size_t i = 0; i < resnames->n_items; ++i) {
        list_t *res_bonds = *((list_t **) dict_get(bonds, resnames->items[i]));

        // order map
        if (maptype == 0) {
            // create an array of pointers to 2D arrays (represented as 1D arrays)
            // each element of this array corresponds to one bond of the residue
            float **bonds_maps = calloc(res_bonds->n_items, sizeof(float *));

            // for each bond, allocate memory for the 2D maps
            for (size_t j = 0; j < res_bonds->n_items; ++j) {
                bonds_maps[j] = calloc(n_tiles, sizeof(float));
            }

            // assign bonds maps to the dictionary
            dict_set(map_dict, resnames->items[i], &bonds_maps, sizeof(float **));
        }

        // samples map
        else if (maptype == 1) {
            size_t **bonds_maps = calloc(res_bonds->n_items, sizeof(size_t *));

            for (size_t j = 0; j < res_bonds->n_items; ++j) {
                bonds_maps[j] = calloc(n_tiles, sizeof(size_t));
            }

            dict_set(map_dict, resnames->items[i], &bonds_maps, sizeof(size_t **));
        }

        else {
            fprintf(stderr, "Incorrect maptype %d (internal error)\n", maptype);
            return NULL;
        }
    }

    return map_dict;
}

void destroy_map(
        dict_t *map, 
        const list_t *resnames, 
        const dict_t *bonds,
        const int maptype)
{
    for (size_t i = 0; i < resnames->n_items; ++i) {
        list_t *res_bonds = *((list_t **) dict_get(bonds, resnames->items[i]));

        if (maptype == 0) {
            float **bonds_maps = *((float ***) dict_get(map, resnames->items[i]));

            for (size_t j = 0; j < res_bonds->n_items; ++j) {
                free(bonds_maps[j]);
            }

            free(bonds_maps);
        }

        else if (maptype == 1) {
            size_t **bonds_maps = *((size_t ***) dict_get(map, resnames->items[i]));

            for (size_t j = 0; j < res_bonds->n_items; ++j) {
                free(bonds_maps[j]);
            }

            free(bonds_maps);
        }

        else {
            fprintf(stderr, "Incorrect maptype %d (internal error)\n", maptype);
            return;
        }
    }

    dict_destroy(map);
}

/* Calculates lipid order maps for one simulation frame. Returns 0 if successful, else returns 1. */
int ordermap_frame(
        atom_selection_t **split,
        const size_t n_residues,
        const dict_t *bonds,
        dict_t *ordermap_upper,
        dict_t *ordermap_lower,
        dict_t *ordermap_full,
        dict_t *samplesmap_upper,
        dict_t *samplesmap_lower,
        dict_t *samplesmap_full,
        const atom_selection_t *phosphates,
        const vec_t membrane_center,
        box_t box,
        const float array_dimx[2],
        const float array_dimy[2],
        const size_t n_cols)
{
    // loop through all the residues
    for (size_t i = 0; i < n_residues; ++i) {
        char *resname = split[i]->atoms[0]->residue_name;

        // get lipid phosphate
        atom_selection_t *phosphate = selection_intersect(split[i], phosphates);
        if (phosphate == NULL || phosphate->n_atoms <= 0) {
            fprintf(stderr, "No phosphate detected for lipid %s (resid %d).\n", resname, split[i]->atoms[0]->residue_number);
            free(phosphate);
            return 1;
        }
        if (phosphate->n_atoms > 1) {
            fprintf(stderr, "Multiple phosphates detected for lipid %s (resid %d).\n", resname, split[i]->atoms[0]->residue_number);
            free(phosphate);
            return 1;
        }

        // assign the lipid into membrane leaflet based on the position of the phosphate
        // 1 -> upper, 0 -> lower
        atom_t *pho = phosphate->atoms[0];
        short classification = 0;
        if (distance1D(pho->position, membrane_center, z, box) > 0) classification = 1;

        free(phosphate);

        dict_t *ordermap_to_use = ordermap_lower;
        dict_t *samplesmap_to_use = samplesmap_lower;
        if (classification) {
            ordermap_to_use = ordermap_upper;
            samplesmap_to_use = samplesmap_upper;
        }

        // get the bonds associated with the residue
        list_t *res_bonds = *((list_t **) dict_get(bonds, resname));
        // for each bond, search for the concerned atoms
        for (size_t b = 0; b < res_bonds->n_items; ++b) {

            char atom_name1[5] = "";
            char atom_name2[5] = "";
            char *bond = list_get(res_bonds, b);
            if (sscanf(bond, "%4s - %4s", atom_name1, atom_name2) != 2) {
                fprintf(stderr, "Could not understand the format of bond %s.\n", bond);
                return 1;
            }

            atom_t *atom1 = NULL;
            atom_t *atom2 = NULL;
            for (size_t j = 0; j < split[i]->n_atoms; ++j) {
                if (strcmp(split[i]->atoms[j]->atom_name, atom_name1) == 0) {
                    if (atom1 != NULL) {
                        fprintf(stderr, "Multiple atoms of the name %s in residue %s.\n", atom_name1, resname);
                        return 1;
                    }

                    atom1 = split[i]->atoms[j];
                    if (atom2 != NULL) break;
                }

                if (strcmp(split[i]->atoms[j]->atom_name, atom_name2) == 0) {
                    if (atom2 != NULL) {
                        fprintf(stderr, "Multiple atoms of the name %s in residue %s.\n", atom_name2, resname);
                        return 1;
                    }

                    atom2 = split[i]->atoms[j];
                    if (atom1 != NULL) break;
                }
            }

            // checking sanity of the result
            if (atom1 == NULL) {
                fprintf(stderr, "Could not find atom %s of bond %s of residue %s.\n", atom_name1, bond, resname);
                return 1;
            }
            
            if (atom2 == NULL) {
                fprintf(stderr, "Could not find atom %s of bond %s of residue %s.\n", atom_name2, bond, resname);
                return 1;
            }

            if (atom1 == atom2) {
                fprintf(stderr, "Bond %s of residue %s is formed by two identical atoms.\n", bond, resname);
                return 1;
            }

            // ignore atoms that are outside of the specified grid
            if (atom1->position[0] < array_dimx[0] || atom1->position[0] > array_dimx[1] ||
                atom1->position[1] < array_dimy[0] || atom1->position[1] > array_dimy[1] ||
                atom2->position[0] < array_dimx[0] || atom2->position[0] > array_dimx[1] ||
                atom2->position[1] < array_dimy[0] || atom2->position[1] > array_dimy[1]) {
                    continue;
            }

            // calculate xy position of the bond
            size_t alloc_atoms = 2;
            atom_selection_t *bond_atoms = selection_create(alloc_atoms);
            selection_add_atom(&bond_atoms, &alloc_atoms, atom1);
            selection_add_atom(&bond_atoms, &alloc_atoms, atom2);
            vec_t bond_center = {0.0};
            center_of_geometry(bond_atoms, bond_center, box);

            free(bond_atoms);
            
            // ignore bonds that are outside of the specified grid
            if (bond_center[0] < array_dimx[0] || bond_center[0] > array_dimx[1] ||
                bond_center[1] < array_dimy[0] || bond_center[1] > array_dimy[1]) {
                    continue;
            }

            // calculate the tile of the array to which the bond belongs
            size_t x_index = coor2index(bond_center[0], array_dimx[0]);
            size_t y_index = coor2index(bond_center[1], array_dimy[0]);
            size_t index = y_index * n_cols + x_index;

            // increase the number of samples at target tile
            size_t **bonds_samples = *((size_t ***) dict_get(samplesmap_to_use, resname));
            ++bonds_samples[b][index];
            bonds_samples = *((size_t ***) dict_get(samplesmap_full, resname));
            ++bonds_samples[b][index];

            // calculate order parameter for the bond
            vec_t vector = {0.0, 0.0, 0.0};
            calc_vector(vector, atom1->position, atom2->position, box);

            float norm2 = (vector[0] * vector[0]) + (vector[1] * vector[1]) + (vector[2] * vector[2]);

            // projection to the z axis
            float projection = vector[2];

            // assign the order parameter to target tile
            float **bonds_order = *((float ***) dict_get(ordermap_to_use, resname));
            bonds_order[b][index] += (projection * projection) / norm2;
            bonds_order = *((float ***) dict_get(ordermap_full, resname));
            bonds_order[b][index] += (projection * projection) / norm2;
        }
    }

    return 0;
}

/* Returns 0 if all maps have been successfully written. Else returns 1. */
int write_map(
        const dict_t *ordermap, 
        const dict_t *samplesmap, 
        const list_t *resnames, 
        const dict_t *bonds,
        const char *output_directory, 
        const char *type,
        const int argc,
        char **argv,
        const size_t n_rows,
        const size_t n_cols,
        const int nan_limit,
        const float minx,
        const float miny)
{
    for (size_t i = 0; i < resnames->n_items; ++i) {
        list_t *res_bonds = *((list_t **) dict_get(bonds, resnames->items[i]));

        float **bonds_maps = *((float ***) dict_get(ordermap, resnames->items[i]));
        size_t **samples_maps = *((size_t ***) dict_get(samplesmap, resnames->items[i]));

        for (size_t j = 0; j < res_bonds->n_items; ++j) {
            // get path to file
            char path_to_file[2048] = "";
            char bond_compact[20] = "";
            memcpy(bond_compact, res_bonds->items[j], strlen(res_bonds->items[j]) + 1);
            strremwhite(bond_compact);
            
            sprintf(path_to_file, "%s/%s_%s_%s.dat", output_directory, type, resnames->items[i], bond_compact);

            // open file
            FILE *output = fopen(path_to_file, "w");
            if (output == NULL) {
                return 1;
            }

            // write header for the file
            fprintf(output, "# Generated with ordermap (C Lipid Order Map Calculator) %s\n", VERSION);
            fprintf(output, "# Command line: ");
            for (int i = 0; i < argc; ++i) {
                fprintf(output, "%s ", argv[i]);
            }
            fprintf(output, "\n");
            fprintf(output, "@ title %s - %s - %s\n", type, resnames->items[i], bond_compact);
            fprintf(output, "@ titlefontsize 18\n");
            fprintf(output, "@ xlabel x coordinate [nm]\n");
            fprintf(output, "@ ylabel y coordinate [nm]\n");
            fprintf(output, "@ zlabel order parameter [arb. u.]\n");
            fprintf(output, "@ zrange -1.0 1.0 0.2\n");
            fprintf(output, "@ grid --\n");
            fprintf(output, "$ type colorbar\n");
            fprintf(output, "$ colormap seismic_r\n");

            // print the map
            for (size_t y = 0; y < n_rows; ++y) {
                for (size_t x = 0; x < n_cols; ++x) {
                    // check that there is enough data for this grid tile
                    if (samples_maps[j][y * n_cols + x] < (size_t) nan_limit) {
                        fprintf(output, "%f %f nan\n", index2coor(x, minx), index2coor(y, miny));
                        continue;
                    }

                    // calculate final order parameter
                    float order = 0.5 * (3.0 * (bonds_maps[j][y * n_cols + x] / samples_maps[j][y * n_cols + x]) - 1.0);

                    // write the final order parameter
                    fprintf(output, "%f %f %f\n", index2coor(x, minx), index2coor(y, miny), order);
                }
            }

            fclose(output);

        }
    }

    return 0;
}


int main(int argc, char **argv)
{
    // get arguments
    char *gro_file = NULL;
    char *xtc_file = NULL;
    char *ndx_file = "index.ndx";
    char *atoms = "Membrane";
    char *itp_file = "martini_v3.0.0_phospholipids_v1.itp";
    char *phosphates = "name PO4";
    float array_dimx[2] = {0.0};
    float array_dimy[2] = {0.0};
    char *output_directory = "order_maps";
    int nan_limit = 30;

    if (get_arguments(argc, argv, &gro_file, &xtc_file, &ndx_file, &itp_file, &atoms, &phosphates, array_dimx, array_dimy, &output_directory, &nan_limit) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    // read gro file
    system_t *system = load_gro(gro_file);
    if (system == NULL) return 1;

    // if array dimensions were not set, get them from gro file
    if (array_dimx[0] == 0 && array_dimx[1] == 0) {
        array_dimx[1] = system->box[0];
    }
    if (array_dimy[0] == 0 && array_dimy[1] == 0) {
        array_dimy[1] = system->box[1];
    }

    // check that the array dimensions don't have nonsensical values
    if (array_dimx[0] >= array_dimx[1] || array_dimy[0] >= array_dimy[1]) {
        fprintf(stderr, "Nonsensical array dimensions.\n");
        free(system);
        return 1;
    }

    // check that no array dimension is negative
    if (array_dimx[0] < 0 || array_dimy[0] < 0) {
        fprintf(stderr, "This program currently does not support negative grid dimensions.\n");
        free(system);
        return 1;
    }

    print_arguments(stdout, gro_file, xtc_file, ndx_file, itp_file, atoms, phosphates, array_dimx, array_dimy, output_directory, nan_limit);

    // check whether the output directory already exists
    DIR *dir = opendir(output_directory);
    if (dir) closedir(dir);
    else if (ENOENT == errno) {
        // if the directory does not exist, create it
        if (mkdir(output_directory, 0750) != 0) {
            fprintf(stderr, "Could not create directory %s.\n", output_directory);
            free(system);
            return 1;
        }
    } else {
        fprintf(stderr, "Could not open or create directory %s.\n", output_directory);
        free(system);
        return 1;
    }

    // try reading ndx file (ignore if this fails)
    dict_t *ndx_groups = read_ndx(ndx_file, system);

    // select all atoms
    atom_selection_t *all = select_system(system);

    // select phosphates
    atom_selection_t *phosphates_select = smart_select(all, phosphates, ndx_groups);
    if (phosphates_select == NULL || phosphates_select->n_atoms == 0) {
        fprintf(stderr, "No phosphates (%s) found.\n", phosphates);

        dict_destroy(ndx_groups);
        free(all);
        free(phosphates_select);
        free(system);
        return 1;
    }

    // select atoms to be analyzed
    select_t *selection = smart_select(all, atoms, ndx_groups);
    if (selection == NULL || selection->n_atoms == 0) {
        fprintf(stderr, "No atoms ('%s') found.\n", atoms);

        dict_destroy(ndx_groups);
        free(selection);
        free(phosphates_select);
        free(all);
        free(system);
        return 1;
    }

    // get residue names to be analyzed
    list_t *resnames = selection_getresnames(selection);

    // get bonds associated with the residues
    dict_t *bonds = get_bonds(itp_file, resnames);
    if (bonds == NULL) {
        fprintf(stderr, "File %s could not be read.\n", itp_file);

        dict_destroy(ndx_groups);
        free(all);
        free(selection);
        free(phosphates_select);
        list_destroy(resnames);
        free(system);
        return 1;
    }

    // check that bond information for all relevant residues have been read
    for (size_t i = 0; i < resnames->n_items; ++i) {
        void *pointer = dict_get(bonds, resnames->items[i]);
        if (pointer == NULL) {
            fprintf(stderr, "Could not find bond information for residue %s.\n", resnames->items[i]);
            
            destroy_bonds(bonds);
            dict_destroy(ndx_groups);
            list_destroy(resnames);
            free(system);
            free(all);
            free(selection);
            free(phosphates_select);
            return 1;
        }
    }

    // split atoms into the individual residues
    atom_selection_t **split = NULL;

    size_t n_residues = selection_splitbyres(selection, &split);
    if (split == NULL || n_residues == 0) {
        fprintf(stderr, "Could not split atoms based on residue number.\n");
        destroy_bonds(bonds);
        list_destroy(resnames);
        
        dict_destroy(ndx_groups);
        free(system);
        free(all);
        free(selection);
        free(phosphates_select);
        return 1;
    }

    // free selection of all atoms
    free(all);

    // TODO: check XTC
    // open xtc file for reading
    XDRFILE *xtc = xdrfile_open(xtc_file, "r");
    if (xtc == NULL) {
        fprintf(stderr, "File %s could not be read as an xtc file.\n", xtc_file);

        destroy_split(split, n_residues);
        destroy_bonds(bonds);
        dict_destroy(ndx_groups);
        list_destroy(resnames);

        free(phosphates_select);
        free(selection);
        free(system);

        return 1;
    }

    if (!validate_xtc(xtc_file, (int) system->n_atoms)) {
        fprintf(stderr, "\nNumber of atoms in %s does not match %s.\n", xtc_file, gro_file);

        destroy_split(split, n_residues);
        destroy_bonds(bonds);
        dict_destroy(ndx_groups);
        list_destroy(resnames);

        free(phosphates_select);
        free(selection);
        free(system);
        xdrfile_close(xtc);
        return 1;
    }

    // prepare dictionaries for lipid order maps
    // ordermaps store information about the lipid order
    size_t n_cols = (size_t) roundf( (array_dimx[1] - array_dimx[0]) * GRID_TILE ) + 1;
    size_t n_rows = (size_t) roundf( (array_dimy[1] - array_dimy[0]) * GRID_TILE ) + 1;
    size_t n_tiles = n_rows * n_cols;

    dict_t *ordermap_upper = create_map(resnames, bonds, n_tiles, 0);
    dict_t *ordermap_lower = create_map(resnames, bonds, n_tiles, 0);
    dict_t *ordermap_full  = create_map(resnames, bonds, n_tiles, 0);

    // samplesmaps store information about the number of samples
    dict_t *samplesmap_upper = create_map(resnames, bonds, n_tiles, 1);
    dict_t *samplesmap_lower = create_map(resnames, bonds, n_tiles, 1);
    dict_t *samplesmap_full  = create_map(resnames, bonds, n_tiles, 1);

    
    int return_code = 0;

    // loop through all the frames in the xtc file
    while (read_xtc_step(xtc, system) == 0) {
        // print info about the progress of reading and writing
        if ((int) system->time % PROGRESS_FREQ == 0) {
            printf("Step: %d. Time: %.0f ps\r", system->step, system->time);
            fflush(stdout);
        }

        // calculate membrane center
        vec_t membrane_center = {0.0};
        center_of_geometry(selection, membrane_center, system->box);

        // calculate lipid order
        if (ordermap_frame(split, n_residues, bonds, ordermap_upper, ordermap_lower, ordermap_full,
            samplesmap_upper, samplesmap_lower, samplesmap_full, phosphates_select,
            membrane_center, system->box, array_dimx, array_dimy, n_cols) != 0) {
                fprintf(stderr, "Could not calculate lipid order for step %d.\n", system->step);
                return_code = 1;
                goto program_end;
        }
    }

    // write maps
    printf("\nWriting maps...\n");
    if (write_map(ordermap_upper, samplesmap_upper, resnames, bonds, output_directory, "upper", 
        argc, argv, n_rows, n_cols, nan_limit, array_dimx[0], array_dimy[0]) != 0) {
            fprintf(stderr, "Could not write maps for upper leaflet.\n");
            return_code = 1;
            goto program_end;
    };

    if (write_map(ordermap_lower, samplesmap_lower, resnames, bonds, output_directory, "lower",
        argc, argv, n_rows, n_cols, nan_limit, array_dimx[0], array_dimy[0]) != 0) {
            fprintf(stderr, "Could not write maps for lower leaflet.\n");
            return_code = 1;
            goto program_end;
    };

    if (write_map(ordermap_full, samplesmap_full, resnames, bonds, output_directory, "full",
        argc, argv, n_rows, n_cols, nan_limit, array_dimx[0], array_dimy[0]) != 0) {
            fprintf(stderr, "Could not write maps for full membrane.\n");
            return_code = 1;
            goto program_end;
    };

    // free memory
    program_end:
    xdrfile_close(xtc);

    destroy_map(ordermap_upper, resnames, bonds, 0);
    destroy_map(ordermap_lower, resnames, bonds, 0);
    destroy_map(ordermap_full,  resnames, bonds, 0);
    
    destroy_map(samplesmap_upper, resnames, bonds, 1);
    destroy_map(samplesmap_lower, resnames, bonds, 1);
    destroy_map(samplesmap_full,  resnames, bonds, 1);

    destroy_split(split, n_residues);
    destroy_bonds(bonds);
    dict_destroy(ndx_groups);
    list_destroy(resnames);
    free(phosphates_select);
    free(selection);
    free(system);

    return return_code;
}