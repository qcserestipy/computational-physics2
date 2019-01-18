program main 
    use cp_mod
    use FT_mod

    character(len=72) :: filename = 'velos.txt'
    type( velocity ), dimension(6000) :: t_velo
 
    call MDReader(filename,t_velo)
end program main