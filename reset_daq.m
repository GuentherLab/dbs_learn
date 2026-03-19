function reset_daq()
    try 
        dio = digitalio('nidaq','Dev1');
        addline(dio, 0:7, 'out');
        data = [0 0 0 0 0 0 0 0];
        putvalue(dio,data);
        fprintf('Dev1 is found, reset done.')
    end

    try 
        dio = digitalio('nidaq','Dev2');
        addline(dio, 0:7, 'out');
        data = [0 0 0 0 0 0 0 0];
            putvalue(dio,data);
        fprintf('Dev2 is found, reset done.')
    end

end