DIAGNOSTIC RUN, not a result row.

The 3.9.1 run in the parent directory dies on the first keyword it does not
know, \energy_term{cmap}, and therefore never reaches the structure.  This run
is the identical input with that one line removed, and nothing else changed,
so that the NEXT thing 3.9.1 cannot do on this system is on the record too.
The \cmap_database_file line is left in place deliberately.
