// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

/*
 * Calculates lipid order parameter for any CG lipids which itp file has been provided.
 * Assumes that the membrane is built in the xy plane with its normal pointing along the z axis.
 */

/* What do the numbers mean?
 * 1    -> perfect alignment with the bilayer normal
 * 0    -> random orientation
 * -0.5 -> anti-alignment
 */

#include <unistd.h>
#include "general.h"

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;

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
        char **phosphates) 
{
    int gro_specified = 0;

    int opt = 0;
    while((opt = getopt(argc, argv, "c:f:n:i:a:p:h")) != -1) {
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
        default:
            return 1;
        }
    }

    if (!gro_specified) {
        fprintf(stderr, "Gro file must always be supplied.\n");
        return 1;
    }
    return 0;
}

void print_usage(const char *program_name)
{
    printf("Usage: %s -c GRO_FILE [OPTION]...\n", program_name);
    printf("\nOPTIONS\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read (optional)\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-i STRING        itp file to read (default: martini_v3.0.0_phospholipids_v1.itp)\n");
    printf("-a STRING        selection of atoms to be used for analysis (default: Membrane)\n");
    printf("-p STRING        identifier of lipid heads (default: name PO4)\n");
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
        const char *phosphates)
{
    fprintf(stream, "\nParameters for Lipid Order calculation:\n");
    fprintf(stream, ">>> gro file:         %s\n", gro_file);
    if (xtc_file == NULL) fprintf(stream, ">>> xtc file:         ----\n");
    else fprintf(stream, ">>> xtc file:         %s\n", xtc_file);
    fprintf(stream, ">>> ndx file:         %s\n", ndx_file);
    fprintf(stream, ">>> itp file:         %s\n", itp_file);
    fprintf(stream, ">>> atoms:            %s\n", atoms);
    fprintf(stream, ">>> phosphates:       %s\n", phosphates);
    fprintf(stream, "\n");
}

void destroy_order(dict_t *order, const list_t *resnames)
{
    for (size_t i = 0; i < resnames->n_items; ++i) {
        float *bonds_order = *((float **) dict_get(order, resnames->items[i]));
        free(bonds_order);
    }
    
    dict_destroy(order);
}

/* Create dictionary for saving order parameters for the individual residue types and bonds.
 * Deallocate using destroy_order().
 */
dict_t *create_order(const list_t *resnames, const dict_t *bonds)
{
    dict_t *order = dict_create();
    for (size_t i = 0; i < resnames->n_items; ++i) {
        list_t *res_bonds = *((list_t **) dict_get(bonds, resnames->items[i]));

        float *bonds_order = calloc(res_bonds->n_items, sizeof(float));
        dict_set(order, resnames->items[i], &bonds_order, sizeof(float *));
    }

    return order;
}

void print_order(const dict_t *order, const list_t *resnames, const dict_t *bonds, const size_t n_frames)
{
    for (size_t i = 0; i < resnames->n_items; ++i) {
        list_t *res_bonds = *((list_t **) dict_get(bonds, resnames->items[i]));
        float *bonds_order = *((float **) dict_get(order, resnames->items[i]));

        printf("> %s <\n", resnames->items[i]);
        for (size_t j = 0; j < res_bonds->n_items; ++j) {
            printf("%s: % .3f\n", res_bonds->items[j], bonds_order[j] / n_frames);
        }
        printf("\n");
    }
}

void add_order(
        dict_t *order_dest, 
        const dict_t *order_src, 
        const list_t *resnames, 
        const dict_t *bonds, 
        const dict_t *n_residues)
{
    for (size_t i = 0; i < resnames->n_items; ++i) {
        list_t *res_bonds = *((list_t **) dict_get(bonds, resnames->items[i]));

        float *src_bonds_order = *((float **) dict_get(order_src, resnames->items[i]));
        float *dest_bonds_order = *((float **) dict_get(order_dest, resnames->items[i]));

        void *raw = dict_get(n_residues, resnames->items[i]);
        // if the entry is not found in the dictionary, that means that no lipids of this type
        // have been dectected and all bond orders are set to 0
        if (raw == NULL) {
            for (size_t j = 0; j < res_bonds->n_items; ++j) {
                dest_bonds_order[j] = 0.0;
            }
            continue;
        }

        // else, calculate bond order
        size_t nres = *( (size_t *) raw);
        for (size_t j = 0; j < res_bonds->n_items; ++j) {
            dest_bonds_order[j] += 0.5 * (3.0 * (src_bonds_order[j] / nres) - 1.0);
        }
    }
}

void increase_n_residues(dict_t *n_residues_dict, const char *resname)
{
    void *raw = dict_get(n_residues_dict, resname);
    
    // if there is no entry for resname in the dictionary, create it
    if (raw == NULL) {
        size_t number = 1;
        dict_set(n_residues_dict, resname, &number, sizeof(size_t));
    // else, increase the number of residues in the entry of a dictionary
    } else {
        *( (size_t *) raw) += 1;
    }
}

static inline void destroy_dictionaries(
        dict_t *order_upper, 
        dict_t *order_lower,
        dict_t *order_full,
        dict_t *n_residues_upper,
        dict_t *n_residues_lower,
        dict_t *n_residues_full,
        const list_t *resnames)
{
    destroy_order(order_upper, resnames);
    destroy_order(order_lower, resnames);
    destroy_order(order_full,  resnames);
    dict_destroy(n_residues_upper);
    dict_destroy(n_residues_lower);
    dict_destroy(n_residues_full);
}

/* Calculates lipid order for one simulation frame. Returns 0 if successful, else returns 1. */
int order_frame(
        atom_selection_t **split,
        const size_t n_residues,
        const dict_t *bonds,
        dict_t *order_upper,
        dict_t *order_lower,
        dict_t *order_full,
        const list_t *resnames,
        const atom_selection_t *phosphates,
        const vec_t membrane_center,
        const box_t box)
{
    // prepare dictionaries for the current lipid order
    dict_t *curr_order_upper = create_order(resnames, bonds);
    dict_t *curr_order_lower = create_order(resnames, bonds);
    dict_t *curr_order_full  = create_order(resnames, bonds);

    // prepare dictionaries for the numbers of residues in each leaflet
    dict_t *n_residues_upper = dict_create();
    dict_t *n_residues_lower = dict_create();
    dict_t *n_residues_full  = dict_create();

    // loop through all the residues
    for (size_t i = 0; i < n_residues; ++i) {
        char *resname = split[i]->atoms[0]->residue_name;

        // get lipid phosphate
        atom_selection_t *phosphate = selection_intersect(split[i], phosphates);
        if (phosphate == NULL || phosphate->n_atoms <= 0) {
            fprintf(stderr, "No phosphate detected for lipid %s (resid %d).\n", resname, split[i]->atoms[0]->residue_number);
            free(phosphate);
            destroy_dictionaries(curr_order_upper, curr_order_lower, curr_order_full, n_residues_upper, n_residues_lower, n_residues_full, resnames);
            return 1;
        }
        if (phosphate->n_atoms > 1) {
            fprintf(stderr, "Multiple phosphates detected for lipid %s (resid %d).\n", resname, split[i]->atoms[0]->residue_number);
            free(phosphate);
            destroy_dictionaries(curr_order_upper, curr_order_lower, curr_order_full, n_residues_upper, n_residues_lower, n_residues_full, resnames);
            return 1;
        }

        // assign the lipid into membrane leaflet based on the position of the phosphate
        // 1 -> upper, 0 -> lower
        atom_t *pho = phosphate->atoms[0];
        short classification = 0;
        if (distance1D(pho->position, membrane_center, z, box) > 0) classification = 1;

        free(phosphate);

        // add lipid to the n_residues dictionaries
        if (classification) increase_n_residues(n_residues_upper, resname);
        else increase_n_residues(n_residues_lower, resname);
        increase_n_residues(n_residues_full, resname);

        // get the bonds associated with the residue
        list_t *res_bonds = *((list_t **) dict_get(bonds, resname));

        // for each bond, search for the concerned atoms
        for (size_t b = 0; b < res_bonds->n_items; ++b) {
            char atom_name1[5] = "";
            char atom_name2[5] = "";
            char *bond = list_get(res_bonds, b);
            if (sscanf(bond, "%4s - %4s", atom_name1, atom_name2) != 2) {
                fprintf(stderr, "Could not understand the format of bond %s.\n", bond);
                destroy_dictionaries(curr_order_upper, curr_order_lower, curr_order_full, n_residues_upper, n_residues_lower, n_residues_full, resnames);
                return 1;
            }

            atom_t *atom1 = NULL;
            atom_t *atom2 = NULL;
            for (size_t j = 0; j < split[i]->n_atoms; ++j) {
                if (strcmp(split[i]->atoms[j]->atom_name, atom_name1) == 0) {
                    if (atom1 != NULL) {
                        fprintf(stderr, "Multiple atoms of the name %s in residue %s.\n", atom_name1, resname);
                        destroy_dictionaries(curr_order_upper, curr_order_lower, curr_order_full, n_residues_upper, n_residues_lower, n_residues_full, resnames);
                        return 1;
                    }

                    atom1 = split[i]->atoms[j];
                    if (atom2 != NULL) break;
                }

                if (strcmp(split[i]->atoms[j]->atom_name, atom_name2) == 0) {
                    if (atom2 != NULL) {
                        fprintf(stderr, "Multiple atoms of the name %s in residue %s.\n", atom_name2, resname);
                        destroy_dictionaries(curr_order_upper, curr_order_lower, curr_order_full, n_residues_upper, n_residues_lower, n_residues_full, resnames);
                        return 1;
                    }

                    atom2 = split[i]->atoms[j];
                    if (atom1 != NULL) break;
                }
            }

            // checking sanity of the result
            if (atom1 == NULL) {
                fprintf(stderr, "Could not find atom %s of bond %s of residue %s.\n", atom_name1, bond, resname);
                destroy_dictionaries(curr_order_upper, curr_order_lower, curr_order_full, n_residues_upper, n_residues_lower, n_residues_full, resnames);
                return 1;
            }
            
            if (atom2 == NULL) {
                fprintf(stderr, "Could not find atom %s of bond %s of residue %s.\n", atom_name2, bond, resname);
                destroy_dictionaries(curr_order_upper, curr_order_lower, curr_order_full, n_residues_upper, n_residues_lower, n_residues_full, resnames);
                return 1;
            }

            if (atom1 == atom2) {
                fprintf(stderr, "Bond %s of residue %s is formed by two identical atoms.\n", bond, resname);
                destroy_dictionaries(curr_order_upper, curr_order_lower, curr_order_full, n_residues_upper, n_residues_lower, n_residues_full, resnames);
                return 1;
            }

            // calculate order parameter for this bond
            vec_t vector = {0.0, 0.0, 0.0};
            calc_vector(vector, atom1->position, atom2->position, box);

            float norm2 = (vector[0] * vector[0]) + (vector[1] * vector[1]) + (vector[2] * vector[2]);

            // projection to the z axis
            float projection = vector[2];

            float *bonds_order = NULL;
            if (classification == 1) bonds_order = *((float **) dict_get(curr_order_upper, resname));
            else bonds_order = *((float **) dict_get(curr_order_lower, resname));

            bonds_order[b] += (projection * projection) / norm2;

            float *bonds_order_full = *((float **) dict_get(curr_order_full, resname));
            bonds_order_full[b] += (projection * projection) / norm2;
        }
    }

    // add current orders to the total orders
    add_order(order_upper, curr_order_upper, resnames, bonds, n_residues_upper);
    add_order(order_lower, curr_order_lower, resnames, bonds, n_residues_lower);
    add_order(order_full, curr_order_full, resnames, bonds, n_residues_full);

    // deallocate memory for the current order and for the n_residues dictionaries
    destroy_dictionaries(curr_order_upper, curr_order_lower, curr_order_full, n_residues_upper, n_residues_lower, n_residues_full, resnames);
    
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

    if (get_arguments(argc, argv, &gro_file, &xtc_file, &ndx_file, &itp_file, &atoms, &phosphates) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    print_arguments(stdout, gro_file, xtc_file, ndx_file, itp_file, atoms, phosphates);

    // read gro file
    system_t *system = load_gro(gro_file);
    if (system == NULL) return 1;

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
    
    // prepare lipid order dictionary
    // contains key->array pairs
    // keys are residue names of all selected lipids
    // arrays contain order parameters for each individual bond (in the same order as in 'bonds')
    dict_t *order_upper = create_order(resnames, bonds);
    dict_t *order_lower = create_order(resnames, bonds);
    dict_t *order_full = create_order(resnames, bonds);

    size_t n_frames = 0;
    // if xtc file has not been provided, analyze just the gro file
    if (xtc_file == NULL) {

        vec_t membrane_center = {0.0};
        center_of_geometry(selection, membrane_center, system->box);

        if (order_frame(split, n_residues, bonds, order_upper, order_lower, order_full, resnames, phosphates_select, membrane_center, system->box) != 0) {
            fprintf(stderr, "Could not calculate lipid order for gro file %s.\n", gro_file);

            destroy_split(split, n_residues);
            destroy_order(order_upper, resnames);
            destroy_order(order_lower, resnames);
            destroy_order(order_full, resnames);
            destroy_bonds(bonds);
            dict_destroy(ndx_groups);
            list_destroy(resnames);

            free(phosphates_select);
            free(selection);
            free(system);

            return 1;
        }

        ++n_frames;
    } else {
        // open xtc file for reading
        XDRFILE *xtc = xdrfile_open(xtc_file, "r");
        if (xtc == NULL) {
            fprintf(stderr, "File %s could not be read as an xtc file.\n", xtc_file);

            destroy_split(split, n_residues);
            destroy_order(order_upper, resnames);
            destroy_order(order_lower, resnames);
            destroy_order(order_full, resnames);
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
            destroy_order(order_upper, resnames);
            destroy_order(order_lower, resnames);
            destroy_order(order_full, resnames);
            destroy_bonds(bonds);
            dict_destroy(ndx_groups);
            list_destroy(resnames);

            free(phosphates_select);
            free(selection);
            free(system);
            xdrfile_close(xtc);
            return 1;
        }

        // loop through all the frames in the xtc file
        while (read_xtc_step(xtc, system) == 0) {
            // print info about the progress of reading and writing
            if ((int) system->time % PROGRESS_FREQ == 0) {
                printf("Step: %d. Time: %.0f ps\r", system->step, system->time);
                fflush(stdout);
            }

            vec_t membrane_center = {0.0};
            center_of_geometry(selection, membrane_center, system->box);

            // calculate lipid order
            if (order_frame(split, n_residues, bonds, order_upper, order_lower, order_full, resnames, phosphates_select, membrane_center, system->box) != 0) {
                fprintf(stderr, "Could not calculate lipid order for step %d.\n", system->step);
                destroy_split(split, n_residues);
                destroy_order(order_upper, resnames);
                destroy_order(order_lower, resnames);
                destroy_order(order_full, resnames);
                destroy_bonds(bonds);
                dict_destroy(ndx_groups);
                list_destroy(resnames);

                free(phosphates_select);
                free(selection);
                free(system);
                xdrfile_close(xtc);

                return 1;
            }

            ++n_frames;
        }

        xdrfile_close(xtc);

    }

    // write the information from the order dictionaries
    printf("\n\nUPPER LEAFLET\n");
    print_order(order_upper, resnames, bonds, n_frames);
    printf("\nLOWER LEAFLET\n");
    print_order(order_lower, resnames, bonds, n_frames);
    printf("\nFULL MEMBRANE\n");
    print_order(order_full, resnames, bonds, n_frames);


    printf("\n");
    destroy_order(order_upper, resnames);
    destroy_order(order_lower, resnames);
    destroy_order(order_full, resnames);
    destroy_split(split, n_residues);
    destroy_bonds(bonds);
    dict_destroy(ndx_groups);
    list_destroy(resnames);

    free(phosphates_select);
    free(selection);
    free(system);
    return 0;
}