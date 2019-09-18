#include "mesh_c.h"

#include <hdf5.h>
#include <iomanip>

// ===========================================================================================
// ===========================================================================================

// ======================================================================================
/// This function prints in hdf format on linear cartesian grids (HEX27->HEX8)
/// the following power quantites: n_cmp=1 -> power; n_cmp=2 power+burnup
void MeshC::print_mesh_med2hdf5(std::string file_name  ///< file name (where to print)
) {  // ===================================================================================
  int n_nodes_mesh_c = (_nx + 1) * (_ny + 1) * (_nz + 1);
  int n_elements_mesh_c = (_nx) * (_ny) * (_nz);

  //   if(n_cmp == 0) std::cout << " DONDRA::print_power_med2hdf5: we have  n_cmp=0 !!!";
  hid_t file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  // file

  hsize_t dimsf[2];
  dimsf[0] = n_nodes_mesh_c;
  dimsf[1] = 1;                                   // dataspace dimension vector dimsf
  hid_t dtsp = H5Screate_simple(2, dimsf, NULL);  // dataspace dtsp
  // coors x ------------------------------------------------------------------
  hid_t dtsetxp = H5Dcreate(
      file, "X1", H5T_NATIVE_DOUBLE, dtsp,
#if HDF5_VERSIONM == 1810
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  hid_t status = H5Dwrite(dtsetxp, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_x[0]);
  H5Dclose(dtsetxp);  // close dataspace
                      // coors y ------------------------------------------------------------------
  dtsetxp = H5Dcreate(
      file, "X2", H5T_NATIVE_DOUBLE, dtsp,
#if HDF5_VERSIONM == 1810
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  status = H5Dwrite(dtsetxp, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_y[0]);
  H5Dclose(dtsetxp);  // close dataspace
                      // coors z ------------------------------------------------------------------
  dtsetxp = H5Dcreate(
      file, "X3", H5T_NATIVE_DOUBLE, dtsp,
#if HDF5_VERSIONM == 1810
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  status = H5Dwrite(dtsetxp, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_z[0]);
  H5Dclose(dtsetxp);  // close dataspace

  // connectivity ------------------------------------------------------------------------
  int* map = new int[n_elements_mesh_c * 8];
  for (int i = 0; i < _nx; i++)
    for (int j = 0; j < _ny; j++)
      for (int k = 0; k < _nz; k++) {
        int idx_elem = i + (_nx)*j + k * (_nx) * (_ny);
        int idx_nodes = i + (_nx + 1) * j + k * (_nx + 1) * (_ny + 1);
        map[idx_elem * 8 + 0] = idx_nodes;
        map[idx_elem * 8 + 1] = idx_nodes + 1;
        map[idx_elem * 8 + 2] = idx_nodes + (_nx + 1) + 1;
        map[idx_elem * 8 + 3] = idx_nodes + (_nx + 1);

        map[idx_elem * 8 + 4] = idx_nodes + (_nx + 1) * (_ny + 1);
        map[idx_elem * 8 + 5] = idx_nodes + 1 + (_nx + 1) * (_ny + 1);
        map[idx_elem * 8 + 6] = idx_nodes + (_nx + 1) + 1 + (_nx + 1) * (_ny + 1);
        map[idx_elem * 8 + 7] = idx_nodes + (_nx + 1) + (_nx + 1) * (_ny + 1);
      }
  dimsf[0] = n_elements_mesh_c * 8;
  dimsf[1] = 1;                             // dataspace dimension vector dimsf
  dtsp = H5Screate_simple(2, dimsf, NULL);  // dataspace dtsp
                                            // connectivity storage
  dtsetxp = H5Dcreate(
      file, "MESH0CONN", H5T_NATIVE_INT, dtsp,
#if HDF5_VERSIONM == 1810
      H5P_DEFAULT, H5P_DEFAULT,
#endif
      H5P_DEFAULT);
  status = H5Dwrite(dtsetxp, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &map[0]);
  H5Dclose(dtsetxp);  // close dataspace
  delete[] map;
  // filling the data vector for hdf5 storage -------------------------------------------
  //   double **pow=new double*[n_cmp];
  //   for(int i=0; i<n_cmp; i++) pow[i]=new double[n_elem*8];
  //
  //   for(int i=0; i<n_cmp; i++) {          // 1-power or 2-power+burnup
  //     int icount=0;
  //     for(int j=0; j<n_elem; j++) {       // reactor elements
  //       for(int isub=0; isub<8; isub++) { // Hex27 -> Hex8
  //         pow[i][icount]=phi.getIJ(j,i);
  //         icount++;
  //       }
  //     }
  //   }
  //   // printing in the hdf5 file  ---------------------------------------------------------
  //
  //   hsize_t dimsf[2]; dimsf[0] =n_elem*8;  dimsf[1] = 1;// dataspace dimension vector dimsf
  //   hid_t dtsp = H5Screate_simple(2, dimsf, NULL);      // dataspace dtsp
  //   for(int i=0; i<n_cmp; i++) {                        // print
  //     hid_t dtsetxp = H5Dcreate(file,name_cmp[i].c_str(),H5T_NATIVE_DOUBLE,dtsp,
  // // #if HDF5_VERSIONM == 1810
  //                               H5P_DEFAULT,H5P_DEFAULT,
  // // #endif
  //                               H5P_DEFAULT);
  //     hid_t status=H5Dwrite(dtsetxp,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,& pow[i][0]);
  //     H5Dclose(dtsetxp); // close dataspace
  //
  //   }
  //   // clean an close ---------------------------------------------------------------------
  //   for(int i=0; i<n_cmp; i++) delete[]pow[i]; delete []pow;
  H5Fclose(file);

  return;
}
// ======================================================================================
