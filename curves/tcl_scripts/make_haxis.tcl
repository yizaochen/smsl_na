proc read_all_pdb_files {start end} {
    mol new $start.pdb type pdb
    for {set index [expr $start + 1]} { $index < [expr $end + 1] } { incr index } {
        mol addfile $index.pdb 0
    }
}