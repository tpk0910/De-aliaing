function bp_comparison(x, tx, a, ta, b, tb, c, tc, folder)
    cd(folder);
    bp_a= load(a).beam_pattern;
    bp_b= load(b).beam_pattern;
    bp_c= load(c).beam_pattern;

    figure;
    plot(x, bp_a); hold on;
    plot(x, bp_b); hold on;
    plot(x, bp_c); hold off;
    axis tight;
    legend(ta, tb, tc);
    title("Beam Pattern Comparison");
    xlabel(tx);
    ylabel("Amplitude (dB)");
end