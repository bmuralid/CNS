
# Development notes

2025-09-16 Tue 11:45 AM


- Currently stuck at `RemakeLevel` subroutine. Not clear what exactly needs to be rebuild and what needs to `fillpatch`ed
- Moving on for now.. need to come back


2025-10-24 Fri 03:59 PM

- Wasted 2-3 hours in just getting the tests running
- Loading the environment variables with `invok` did not work
- This resulted in linking with the wrong AMReX library
- Need a foolproof way to set the correct CMAKE root dirs for the various libraries
