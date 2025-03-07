function f_PSD_Average(v_Data, srt, srate, GF, h, Ch)

    srate = srate / srt;
    v_Data = downsample(v_Data, srt);

    [pxx, f] = pwelch(v_Data, [], [], [], srate);

    if GF == 1
        pxx_av = pxx;
    else
        pxx_av = mean(pxx');
    end

    figure;
    plot(f, pxx_av)
    xlabel('Frequency (Hz)')
    ylabel('PSD Average');
    title(h.recChNames(Ch));

end
