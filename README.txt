%%% DendroScan -- Software for multivariate statistical particle shape analyses 
%   after Dürig et al. (2012, 2020)
%
%  Copyright (C) 2020 by Tobias Dürig and Louise Steffensen Schmidt
%
%  Authors: T. Dürig, Institute of Earth Sciences, University of Iceland,
%          Reykjavík, Iceland
%           L. S. Schmidt, Department of Geosciences, University of Oslo,
%           Oslo, Norway
%
%  This program is free software licensed under the GPL (>=v3).
%  Read the GPS-3.TXT file that comes with this software for details.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This program is free software; you can redistribute it and/or modify it
%  under the terms of the GNU General Public License as published by the
%  Free Software Foundation; either version 3 of the License, or (at your
%  option) any later version.
%  
%  Parts of DendroScan are not copyright by the DendroScan development team.
%  The original authors hold the copyrights and you have to abide to their
%  licensing terms where noted. See the headers of the respective .m scripts
%  for details.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License (GPL) for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, see <http://www.gnu.org/licenses/>.
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input is in the form of flat ASCII files in '.csv' format, containing 
%       a list of morphometric parameters. DendroScan was designed to
%       digest tables produced by the software PARTIcle Shape ANalyzer 
%       (PARTISAN)by Dürig et al. (2018).
%      
%  Output files according to selection in the GUI.
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  References:
%   Dürig et al. (2012): doi.org/10.1007/s00445-011-0562-0
%   Dürig et al. (2018): doi.org/10.4401/ag-7865
%   Dürig et al. (2020): in review
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you wish to contribute to the development of DendroScan,
% or to report bugs or other problems with the software, please contact me 
% per email (tobi@hi.is).


System requirements
===================
DendroScan was developed and tested with Matlab R2019b on both Linux (debian) and Windows operation systems. DendroScan can be run with Matlab in default installation and does not require any additional Matlab libraries.

Installation and start
======================
Unzip the software package in a folder of your choice. The program is started by opening “main.m” with Matlab and executing it. 



To install DendroScan, unzip all contents of "DendroScan2_10.zip" in your work directory.
To start the program, run "main.m".


