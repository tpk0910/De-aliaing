cfg = struct;
cfg.cpwc_dir = "C:\Users\mark0\OneDrive\Desktop\菸酒生\De-aliaing\CPWC_dual";
cfg.cpwc_entry = "cpwc_dealiasing";
cfg.pw_dir = "C:\Users\mark0\OneDrive\Desktop\菸酒生\De_aliasing_MingChi\1.5D";
cfg.pw_entry = "PW_dealiasing";

cfg.save_dir = "C:\Users\mark0\OneDrive\Desktop\菸酒生\De-aliaing\Compare";
cfg.showFigures = true;
cfg.alpha = 0.5;

T = compare_cpwc_vs_pw(cfg);

evalin('base','clear;');  % optional: clean base workspace

