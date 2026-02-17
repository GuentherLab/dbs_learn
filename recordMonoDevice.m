function recordMonoDevice(deviceName, outFile)
    fs = 48000;
    frameSize = 4096;

    d = audioDeviceReader( ...
        'Device', deviceName, ...
        'SampleRate', fs, ...
        'NumChannels', 2, ...
        'SamplesPerFrame', frameSize);

    w = dsp.AudioFileWriter(outFile, 'SampleRate', fs);

    fprintf('Recording from %s -> %s\n', deviceName, outFile);

    while true
        x = d();          % Nx2
        mono = x(:,1);    % take left channel only
        w(mono);
    end
end
