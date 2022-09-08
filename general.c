// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include "general.h"

void destroy_bonds(dict_t *bonds)
{
    char **keys = NULL;
    size_t n_keys = dict_keys(bonds, &keys);

    for (size_t i = 0; i < n_keys; ++i) {
        list_t *list = *((list_t **) dict_get(bonds, keys[i]));
        list_destroy(list);
    }

    free(keys);
    dict_destroy(bonds);
}

void destroy_split(atom_selection_t **split, size_t n_residues)
{
    for (size_t i = 0; i < n_residues; ++i) {
        free(split[i]);
    }
    
    free(split);
}

dict_t *get_bonds(char *itp_file, list_t *resnames)
{
    // read itp file and obtain information about bonds
    FILE *itp = fopen(itp_file, "r");
    if (itp == NULL) {
        return NULL;
    }

    dict_t *bonds_dict = dict_create();

    char line[1024] = "";
    int molname = 0;
    int atoms = 0;
    int bonds = 0;
    int block = -1;
    list_t *atom_names = list_create();
    list_t *bond_list = list_create();
    while (fgets(line, 1024, itp) != NULL) {
        // remove gromacs comments
        char *chrp = NULL;
        if ((chrp = strchr(line, ';')) != NULL) {
            *chrp = 0;
        }
        // remove gromacs preprocessor marks
        if ((chrp = strchr(line, '#')) != NULL) {
            *chrp = 0;
        }

        // ignore empty lines
        char **split = NULL;
        char line_copy[1024] = "";
        strncpy(line_copy, line, 1024);
        line_copy[1023] = 0;
        if (strsplit(line_copy, &split, " \t\n") <= 0) continue;
        
        // read residue name on the next non-empty line
        if (strstr(line, "moleculetype") != NULL) {
            molname = 1;
        }

        // check if residue name matches any of the residue names requested for the analysis
        else if (molname) { 
            for (size_t i = 0; i < resnames->n_items; ++i) {
                if (strcmp(split[0], resnames->items[i]) == 0) {
                    block = i;
                    break;
                }
            }

            molname = 0;
        } 

        else if (block != -1) {
            if (strstr(line, "atoms") != NULL) {
                atoms = 1;
            }

            else if (atoms) {
                if (strstr(line, "bonds") != NULL) {
                    atoms = 0;
                    bonds = 1;
                    free(split);
                    continue;
                }

                list_append(&atom_names, split[4]);
            }

            else if (bonds) {
                if (strstr(line, "[") != NULL) {
                    dict_set(bonds_dict, resnames->items[block], &bond_list, sizeof(list_t *));
                    
                    // create new lists for any other residues
                    list_destroy(atom_names);
                    atom_names = list_create();
                    bond_list = list_create();
                    
                    bonds = 0;
                    block = -1;
                    free(split);
                    continue;
                }

                int a1 = 0;
                int a2 = 0;
                if (sscanf(split[0], "%d", &a1) != 1) {
                    free(split);
                    list_destroy(atom_names);
                    list_destroy(bond_list);
                    destroy_bonds(bonds_dict);
                    fclose(itp);
                    return NULL;
                }

                if (sscanf(split[1], "%d", &a2) != 1) {
                    free(split);
                    list_destroy(atom_names);
                    list_destroy(bond_list);
                    destroy_bonds(bonds_dict);
                    fclose(itp);
                    return NULL;
                }

                char format[20] = "";
                sprintf(format, "%s - %s", atom_names->items[a1 - 1], atom_names->items[a2 - 1]);

                list_append(&bond_list, format);
            }

        }

        free(split);
    }

    list_destroy(bond_list);

    fclose(itp);
    list_destroy(atom_names);

    return bonds_dict;
}