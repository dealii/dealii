#pragma once

#include <deal.II/base/exception_macros.h>

#include <boost/algorithm/string/join.hpp>

#include <netcdf.h>

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

void
dump_attribute(LogStream &log, int ncid, int varid, int attid)
{
  int  ncerr;
  char varname[NC_MAX_NAME];
  char attname[NC_MAX_NAME];
  if (varid == NC_GLOBAL)
    {
      strncpy(varname, "global", NC_MAX_NAME);
    }
  else
    {
      ncerr = nc_inq_varname(ncid, varid, varname);
      AssertThrowNC(ncerr);
    }
  ncerr = nc_inq_attname(ncid, varid, attid, attname);
  AssertThrowNC(ncerr);
  nc_type atttype;
  ncerr = nc_inq_atttype(ncid, varid, attname, &atttype);
  AssertThrowNC(ncerr);
  size_t attlen;
  ncerr = nc_inq_attlen(ncid, varid, attname, &attlen);
  AssertThrowNC(ncerr);
  if (atttype == NC_CHAR)
    {
      std::string attval(attlen, char(0));
      ncerr = nc_get_att_text(ncid, varid, attname, attval.data());
      AssertThrowNC(ncerr);
      log << varname << ":" << attname << " = " << attval << std::endl;
    }
  else if (attlen == 1 && atttype == NC_INT)
    {
      int attval;
      ncerr = nc_get_att_int(ncid, varid, attname, &attval);
      AssertThrowNC(ncerr);
      log << varname << ":" << attname << " = " << attval << std::endl;
    }
  else if (attlen == 1 && atttype == NC_DOUBLE)
    {
      double attval;
      ncerr = nc_get_att_double(ncid, varid, attname, &attval);
      AssertThrowNC(ncerr);
      log << varname << ":" << attname << " = " << attval << std::endl;
    }
  else
    {
      log << varname << ":" << attname << " = "
          << "Unexpected attribute type" << std::endl;
    }
}

void
dump_variable(LogStream &log, int ncid, int varid)
{
  char             varname[NC_MAX_NAME];
  std::vector<int> dims(NC_MAX_DIMS);
  int              ndims;
  int              natts;
  int              ncerr;
  ncerr =
    nc_inq_var(ncid, varid, varname, nullptr, &ndims, dims.data(), &natts);
  AssertThrowNC(ncerr);
  dims.resize(ndims);
  std::vector<std::string> dimlens(ndims);
  std::transform(dims.begin(), dims.end(), dimlens.begin(), [&](auto dim) {
    size_t len;
    ncerr = nc_inq_dimlen(ncid, dim, &len);
    AssertThrowNC(ncerr);
    return std::to_string(len);
  });
  log << varname << "(" << boost::join(dimlens, ", ") << ")" << std::endl;
  for (int i = 0; i < natts; ++i)
    {
      dump_attribute(log, ncid, varid, i);
    }
}

void
dump_nc_file(LogStream &log, const char *ncfile)
{
  AssertThrow(std::filesystem::is_regular_file(ncfile),
              ExcIO("Output file not found."));
  int ncerr;
  int ncid;
  ncerr = nc_open(ncfile, 0, &ncid);
  AssertThrowNC(ncerr);
  int ndims, nvars, natts;
  ncerr = nc_inq(ncid, nullptr, &nvars, &natts, nullptr);
  AssertThrowNC(ncerr);
  for (int i = 0; i < natts; ++i)
    {
      dump_attribute(log, ncid, NC_GLOBAL, i);
    }
  std::vector<int> vars(nvars);
  ncerr = nc_inq_varids(ncid, nullptr, vars.data());
  AssertThrowNC(ncerr);
  for (auto var : vars)
    {
      dump_variable(log, ncid, var);
    }
  nc_close(ncid);
}
