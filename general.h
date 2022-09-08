// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#ifndef GENERAL_H
#define GENERAL_H

#include <groan.h>

/*! @brief Destroys bonds dictionary. */
void destroy_bonds(dict_t *bonds);

/*! @brief Destroys an array of atom selection structures. */
void destroy_split(atom_selection_t **split, size_t n_residues);

/*! @brief Reads itp file searching for bonds of relevant residues. */
dict_t *get_bonds(char *itp_file, list_t *resnames);

#endif /* GENERAL_H */