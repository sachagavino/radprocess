/******************************************************************************************
*   This file is part of PyMSES.                                                          *
*                                                                                         *
*   PyMSES is free software: you can redistribute it and/or modify                        *
*   it under the terms of the GNU General Public License as published by                  *
*   the Free Software Foundation, either version 3 of the License, or                     *
*   (at your option) any later version.                                                   *
*                                                                                         *
*   PyMSES is distributed in the hope that it will be useful,                             *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of                        *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                         *
*   GNU General Public License for more details.                                          *
*                                                                                         *
*   You should have received a copy of the GNU General Public License                     *
*   along with PyMSES.  If not, see <http://www.gnu.org/licenses/>.                       *
******************************************************************************************/
#include <stdio.h>
#include "read_amr.h"

int main() {


	char fname[1024];
	int icpu;
	for (icpu=0; icpu<512; icpu++) {
		RAMSES_AmrHeader hdr;
		RAMSES_AmrStruct amr;
		sprintf(fname, "/data/milkyway/smooth/output_00700/amr_00700.out%05d", icpu+1);
		printf("Reading AMR %s\n", fname);
		RAMSES_ReadAmr(&hdr, &amr, fname, 99, 0);
		RAMSES_AmrHeader_Free(&hdr);
		RAMSES_AmrStruct_Free(&amr);
	}

	return 0;
}
