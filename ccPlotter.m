% correlation coefficients plotter

clearvars;

% environment variables
if ismac
    profile='HOME';
else
    profile='USERPROFILE';
end
user_path=getenv(profile);
desktop_path=[user_path,filesep,'Desktop'];
project_path=[desktop_path,filesep,'DyNaVoiceR_Files'];
vers=['Release R' version('-release')]; skip=0;

% powerodft parameters and variables
N       = 1024; N2 = N/2; N4 = N/4;
win     = sin(pi/N*((0:N-1) + 0.5));
FS      = 22050;
data   = zeros(1,N);	% input data
idata   = zeros(1,N);	% input (int) data
fdata   = zeros(1,N);	% (float complex) data
% complex vectors for the non-optimal computation of the ODFT, and IODFT transforms
direxp = exp(-1i*pi*[0:(N-1)]/N);
invexp = exp( 1i*pi*[0:(N-1)]/N);

% input audio file via UI (original .wav file)
[filename, path, ~] = uigetfile([project_path,filesep,'*.wav'], 'WAV Files (*.wav)', 'MultiSelect','on');

% input template files via UI (.mat file)
[t_filename, t_path, ~] = uigetfile([project_path,filesep,'*.mat'], 'MAT Files (*.mat)', 'MultiSelect','on');

n_vowels=length(t_filename);
% to do -> load t_filename -> load model 
figure;
hold on;
for vowel_index=1:n_vowels
    load([t_path,t_filename{vowel_index}]);
    v_model(vowel_index,:)=vowel_model;
    plot(vowel_model);
end
hold off;


for item=1:length(filename)

%get powerodft per frame -> frame model
[audio_original, fs] = audioread([path,filename{item}]);
audio_original=audio_original.';

%zero padding
total_frames=floor(length(audio_original)/N2)+2;
audio_padded(N2+1:N2+length(audio_original))=audio_original;
audio_padded((total_frames+1)*N2)=0;
 
tmpdata=[];
%get powerodft -> powerodft matrix
for step=1:total_frames
     
    data=audio_padded(1+(step-1)*N2:N+(step-1)*N2);
    %get powerodft(step)
    idata=data.*win;
    fdata = idata.*direxp; % this is sub-optimal ODFT computation
    odft=fft(fdata);       % this is sub-optimal ODFT computation
    %phaseodft = angle(odft(1:N2));
    frame_powerodft(step,:) = 1E-6+abs(odft).^2;
    
    frame_model(step,:)=CoSMo(frame_powerodft(step,:));
       
end

% % plot models
% figure;
% hold on
% for step=1:total_frames
%      
%     plot(frame_model(step,:));
%        
% end
% hold off


% to do -> calc cc's 1 x n_models
for step=1:total_frames
    
    for vowel_index=1:n_vowels
        cc(vowel_index,step)=corr(frame_model(step,:).',v_model(vowel_index,:).','type','pearson');
        if step==1
            ccave(vowel_index,step)=cc(vowel_index,step);
        else
            ccave(vowel_index,step)=0.5*cc(vowel_index,step)+0.5*cc(vowel_index,step-1);
        end
    end
       
end

for step=1:total_frames % top results filter
    
    sorted=sort(ccave(2:n_vowels,step),'descend');
    
      for vowel_index=2:n_vowels

            if cc(vowel_index,step)<cc(1,step)
               cc(vowel_index,step)=0; 
            end
            if ( ccave(vowel_index,step)<ccave(1,step) || ccave(vowel_index,step)<sorted(3) )
               ccave(vowel_index,step)=0; 
            end

      end
end

vowel_name_display='siêéáâóôueãêîõû';

% to do -> plot cc's
figure;
hold on;
title([filename{item} ' (raw)']);
line([1 length(cc(1,:))], [0.8 0.8],'Color','black','LineStyle',':');
xlim([0 1+length(cc(1,:))])
for vowel_index=1:n_vowels
    plot(cc(vowel_index,:),'DisplayName',vowel_name_display(vowel_index));
end
hold off;

figure;
hold on;
title([filename{item} ' (averaged)']);
line([1 length(cc(1,:))], [0.8 0.8],'Color','black','LineStyle',':');
xlim([0 1+length(cc(1,:))])
for vowel_index=1:n_vowels
    plot(ccave(vowel_index,:),'DisplayName',vowel_name_display(vowel_index));
end
hold off;

end

disp('done!');