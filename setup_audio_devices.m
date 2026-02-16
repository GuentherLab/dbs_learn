%%%% modify this function based on the task computer and audio devices you are using
% some computers have specifications for input devices; however we might not end up doing audio recording in the caller script

function aud = setup_audio_devices()

aud = struct; 

computername = getenv('COMPUTERNAME'); % might not work on non-windows machines
    computername = deblank(computername); 
auddevs = audiodevinfo; 
    devs_in = {auddevs.input.Name};
    devs_out = {auddevs.output.Name};

% % % %  % default device names do not have to be full device names, just need to be included in a single device name
% % % % if any(contains(devs_out,'Focusrite') )  % Full experimental setup with Focusrite
% % % %         aud.device_in = 'Focusrite';
% % % %         aud.device_out = 'Focusrite';
% % % % elseif any(contains(devs_out,'Speakers (Logitech G432 Gaming Headset)') ) % this name might instead be 'Speakers (G432 Gaming Headset)' if GHub is not installed
% % % %     % if G432 is available, use it for input and output...
% % % %     % .... if using Eprom instead, use 3.5mm jack  
% % % %     aud.device_in = 'Microphone (Logitech G432 Gaming Headset)';
% % % %     aud.device_out = 'Speakers (Logitech G432 Gaming Headset)'; 
% % % % else
    switch computername
        case '677-GUE-WL-0010'  % AM work laptop - Thinkpad X1
            if any(contains(devs_out,'Headphones (WF-C500)') ) % if using bluetooth headphones

                aud.device_in = 'Default';
                aud.device_out = 'Headphones (WF-C500)'; 
                    % aud.device_out = 'Headset (WF-C500)'; 
            else % Thinkpad X1 without headphones
                aud.device_in = 'Microphone'; 
                aud.device_out = 'Realtek'; 
                    % aud.device_out = 'ARZOPA'; % portable screen speakers
            end
        case 'AMSMEIER' % AM strix laptop
            if any(contains(devs_out,'Speakers (Realtek(R) Audio) (Windows DirectSound)')) 
                aud.device_out = 'Speakers (Realtek(R) Audio)';
                aud.device_in = 'Microphone Array (Realtek(R) Audio)'; 
            end

        % % % case {'MSI','amsmeier'} % AM laptop
        % % %     if any(contains(devs_in, 'EEPROM')) % wired headset - usb input
        % % %         aud.device_in = 'EEPROM'; 
        % % %         aud.device_out = 'Realtek HD Audio 2nd'; % if this usb input is being used, output will be via (backmost) 3.5mm jack
        % % %     elseif any(contains(devs_out,'MP43250') )   % if using bluetooth headphones
        % % %         aud.device_ing = 'MP43250'; 
        % % %         aud.device_out = 'MP43250'; 
        % % %     else
        % % %         aud.device_in = 'Microphone (Realtek(R) Audio)';    % onboard mic
        % % %         % aud.device_out = 'Speakers (Realtek(R) Audio)'; % use default output - onboard speakers
        % % %         aud.device_out = 'Headphones (Realtek(R) Audio)'; 
        % % %     end
        
        otherwise 
            error('unknown computer; please add preferred devices to "audio device section" of flvoice_run.m')
    end
% % % end


