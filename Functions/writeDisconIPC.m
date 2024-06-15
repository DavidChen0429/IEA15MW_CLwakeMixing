function writeDisconIPC(IPCenable,Pulse,HelixCM2,HelixCM3,Frequency,PitchOffset)
fid = fopen('discon.in', 'wt');
fprintf(fid,'%.4f	= Uservar1: IPC enable(1) disable (0)\n', IPCenable);
fprintf(fid,'%.4f	= Uservar2: Pulse (amplitude in rad/s and div. by 2 0.0698)\n', Pulse);
fprintf(fid,'%.4f	= Uservar3: Helix CM2 (amplitude in rad and div. by 2 0.0698)\n', HelixCM2);
fprintf(fid,'%.4f	= Uservar4: Helix CM3 (amplitude in rad and div. by 2 0.0698)\n', HelixCM3);
fprintf(fid,'%.4f	= Uservar5: Frequency (rad/s)\n', Frequency);
fprintf(fid,'%.4f	= Uservar6: Pitch Offset (rad)\n', PitchOffset);
fclose(fid);
end
