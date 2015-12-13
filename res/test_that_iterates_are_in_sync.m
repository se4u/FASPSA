cd ~/Dropbox/paper/faspsa/res/
load sso_project_intermediate.mat
s2 = results_struct.FeedbackAdaptive2SPSA_10_2_120000_loss_sequence;
s1 = (results_struct.EfficientFeedbackAdaptive2SPSA_10_2_120000_loss_sequence);
hold off; plot(s1(1:end), 'b'); hold on; plot(s2(1:end), 'r');
pause;
s1 = (results_struct.EfficientFeedbackAdaptive2SPSA_10_1_120000_loss_sequence);
s2 = (results_struct.FeedbackAdaptive2SPSA_10_1_120000_loss_sequence);
hold off; plot(s1(1:end), 'b'); hold on; plot(s2(1:end), 'r');
