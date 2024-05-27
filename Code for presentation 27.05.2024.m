
%%%%% CSI DOWNLINK REPORT DEMO 5G %%%%
% 27.05.2024

% Kluczowe Skróty w modelu CSI Reporting

% CSI Reporting (Channel State Information Reporting) - Proces raportowania CSI polega na zbieraniu informacji o stanie kanału przez urządzenie użytkownika (UE) i 
% przekazywaniu ich do stacji bazowej (gNB) w formie opinii zwrotnej. 
% Wskazuje to, jak dobrze sygnał jest odbierany i interpretowany przez odbiornik. Na podstawie tych raportów gNB może dostosowywać parametry transmisji, aby zoptymalizować jakość i wydajność komunikacji.
% Wszystkie te parametry - RI, PMI, CQI - są zbierane przez UE i przekazywane gNB jako część procesu raportowania CSI. Dzięki tym informacjom gNB może inteligentnie zarządzać transmisją danych, aby zapewnić optymalną wydajność sieci.



%  RI (Rank Indicator) - Wskaźnik rangi określa liczbę możliwych warstw (strumieni danych) dla transmisji w kierunku dolnym, co odpowiada maksymalnej liczbie nieskorelowanych ścieżek, 
% które mogą być wykorzystane do transmisji. Jest to istotne dla gNB, aby optymalnie dostosować liczbę warstw transmisji do aktualnych warunków kanału.

% PMI (Precoding Matrix Indicator) - Wskaźnik macierzy precodowania definiuje odpowiednią macierz przetwarzania sygnału dla transmisji w kierunku dolnym. 
% Wybór odpowiedniej macierzy precodowania jest kluczowy dla zapewnienia optymalnej jakości transmisji, szczególnie w zmiennych warunkach kanału.

% CQI (Channel Quality Indicator) - Wskaźnik jakości kanału informuje o jakości kanału komunikacyjnego na podstawie pomiarów parametrów fizycznych, takich jak SNR 
% (Stosunek sygnału do szumu) czy SINR (Stosunek sygnału do interferencji plus szum). 
% Jest to istotne dla gNB do dynamicznego dostosowywania parametrów transmisji, takich jak modulacja i szybkość kodowania.



% Liczba 10 ms ramek
nFrames = 5; % Określenie liczby ramek na 10 ms

% SNR w dB
SNRdB = 10; % Ustawienie wartości stosunku sygnału do szumu w decybelach

% Konfiguracja nośnika
carrier = nrCarrierConfig; % Inicjalizacja konfiguracji nośnika
carrier.SubcarrierSpacing = 15; % Ustawienie odstępu między nośnymi
carrier.NSizeGrid = 52; % Określenie rozmiaru siatki nośnej

% Parametry BWP
NStartBWP = 0; % Określenie początkowego numeru pasma bazowego
NSizeBWP = 52; % Określenie rozmiaru pasma bazowego

% Konfiguracja CSI-RS

csirs = nrCSIRSConfig; % Inicjalizacja konfiguracji CSI-RS CSI-RS (Channel State Information Reference Signal) 
% to specjalne sygnały referencyjne używane w sieciach komórkowych (np. 5G) do oceny stanu kanału. 

csirs.CSIRSType = {'nzp','nzp','nzp'}; % Określenie typu sygnału referencyjnego CSI-RS
csirs.RowNumber = [4 4 4]; % Określenie liczby wierszy
csirs.NumRB = 52; % Określenie liczby pasm bazowych
csirs.RBOffset = 0; % Określenie przesunięcia pasma bazowego
csirs.CSIRSPeriod = [4 0]; % Określenie okresu sygnału referencyjnego CSI-RS
csirs.SymbolLocations = {0, 0, 0}; % Określenie lokalizacji symboli
csirs.SubcarrierLocations = {0, 4, 8}; % Określenie lokalizacji nośnych
csirs.Density = {'one','one','one'}; % Określenie gęstości sygnału referencyjnego CSI-RS

% Liczba anten nadawczych i odbiorczych, musi byc tyle samo ich 
nTxAnts = csirs.NumCSIRSPorts(1); % Określenie liczby anten nadawczych
nRxAnts = 4; % Określenie liczby anten odbiorczych

% Walidacja konfiguracji CSI-RS
validateCSIRSConfig(carrier,csirs,nTxAnts); % Sprawdzenie poprawności konfiguracji CSI-RS



            % CQITable (Tabela CQI):
            
            % reportConfig.CQITable = 'table1';
            % Wybór tabeli CQI (Channel Quality Indicator) używanej do mapowania jakości kanału na odpowiednie wskaźniki CQI. 
            % Przykładowo, tabela CQI określa, jakie wartości CQI są używane do określenia modulacji i kodowania (MCS, Modulation and Coding Scheme).
            
            % CodebookType (Typ tablicy kodów):
            % reportConfig.CodebookType = 'Type1SinglePanel';
            % Określa typ tablicy kodów, który jest używany do precodingu w systemie MIMO (Multiple Input Multiple Output). 
            % Typ 1 oznacza podstawowy typ tablicy kodów, a 'SinglePanel' wskazuje na jednopanelową strukturę.
            % 

% Konfiguracja raportowania CSI
reportConfig.NStartBWP = NStartBWP; % Ustawienie początkowego numeru pasma bazowego w raporcie
reportConfig.NSizeBWP = NSizeBWP; % Ustawienie rozmiaru pasma bazowego w raporcie
reportConfig.CQITable = 'table1'; % Określenie tabeli CQI
reportConfig.CodebookType = 'Type1SinglePanel'; % Określenie typu tablicy kodów
reportConfig.PanelDimensions = [2 1]; % Określenie wymiarów panelu
reportConfig.CQIMode = 'Subband'; % Określenie trybu CQI
reportConfig.PMIMode = 'Subband'; % Określenie trybu PMI
reportConfig.SubbandSize = 4; % Określenie rozmiaru podpasma
reportConfig.PRGSize = []; % Określenie rozmiaru PRG ,PRG (Physical Resource Group): Grupa bloków PRB, używana do alokacji zasobów w ramach 5G NR (New Radio).
reportConfig.CodebookMode = 1; % Określenie trybu tablicy kodów
reportConfig.CodebookSubsetRestriction = []; % Określenie ograniczenia podzbioru tablicy kodów
reportConfig.i2Restriction = []; % Określenie ograniczenia indeksu i2
reportConfig.RIRestriction = []; % Określenie ograniczenia RI
reportConfig.NumberOfBeams = 2; % Określenie liczby wiązek
reportConfig.SubbandAmplitude = false; % Określenie amplitudy podpasma
reportConfig.PhaseAlphabetSize = 4; % Określenie rozmiaru alfabetu fazowego
reportConfig.ParameterCombination = 2; % Określenie kombinacji parametrów
reportConfig.NumberOfPMISubbandsPerCQISubband = 2; % Określenie liczby podpasm PMI na podpasma CQI

% Konfiguracja kanału
channel = nrTDLChannel; % Inicjalizacja kanału
channel.NumTransmitAntennas = nTxAnts; % Określenie liczby anten nadawczych
channel.NumReceiveAntennas = nRxAnts; % Określenie liczby anten odbiorczych
channel.DelayProfile = 'TDL-C'; % Określenie profilu opóźnienia
channel.MaximumDopplerShift = 50; % Określenie maksymalnego przesunięcia Dopplera
channel.DelaySpread = 300e-9; % Określenie rozpraszania opóźnienia

% Informacje o OFDM
waveformInfo = nrOFDMInfo(carrier); % Uzyskanie informacji o OFDM
channel.SampleRate = waveformInfo.SampleRate; % Ustawienie częstotliwości próbkowania kanału

% Informacje o kanale
chInfo = info(channel); % Uzyskanie informacji o kanale
maxChDelay = chInfo.MaximumChannelDelay; % Uzyskanie maksymalnego opóźnienia kanału

% Całkowita liczba slotów
totSlots = nFrames * carrier.SlotsPerFrame; % Obliczenie całkowitej liczby slotów

% Inicjalizacja zmiennych do przechowywania wyników
cqiPracticalPerSlot = []; % Inicjalizacja macierzy do przechowywania praktycznych wartości CQI
subbandCQIPractical = []; % Inicjalizacja macierzy do przechowywania praktycznych wartości podpasma CQI
pmiPracticalPerSlot = struct('i1', [], 'i2', []); % Inicjalizacja struktury do przechowywania praktycznych wartości PMI
SINRPerSubbandPerCWPractical = []; % Inicjalizacja macierzy do przechowywania praktycznych wartości SINR
cqiPerfectPerSlot = []; % Inicjalizacja macierzy do przechowywania idealnych wartości CQI
subbandCQIPerfect = []; % Inicjalizacja macierzy do przechowywania idealnych wartości podpasma CQI
pmiPerfectPerSlot = struct('i1', [], 'i2', []); % Inicjalizacja struktury do przechowywania idealnych wartości PMI
SINRPerSubbandPerCWPerfect = []; % Inicjalizacja macierzy do przechowywania idealnych wartości SINR
riPracticalPerSlot = []; % Inicjalizacja wektora do przechowywania praktycznych wartości RI
riPerfectPerSlot = []; % Inicjalizacja wektora do przechowywania idealnych wartości RI

% Liczba portów CSI-RS
csirsPorts = csirs.NumCSIRSPorts(1); % Pobranie liczby portów CSI-RS

% Długości CDM (Code Division Multiplexing )dla skonfigurowanych zasobów CSI-RS
cdmLengths = getCDMLengths(csirs); % Obliczenie długości CDM dla zasobów CSI-RS

% Inicjalizacja przesunięcia czasowego
offsetPractical = 0; % Ustawienie początkowego przesunięcia czasowego

% NZP-CSI-RS (Non-Zero Power Channel State Information Reference Signal)
% NZP-CSI-RS są sygnałami referencyjnymi, które mają niezerową moc i są używane do uzyskiwania dokładnych pomiarów stanu kanału.

totSlotsBinaryVec = zeros(1, totSlots); % Inicjalizacja wektora binarnego informującego o obecności NZP-CSI-RS w slotach

% Ustawienie stanu RNG dla powtarzalności
rng('default'); % Ustawienie stanu generatora liczb losowych na to samo dla powtarzalnoscci 

% Pętla po wszystkich slotach
for nslot = 0:totSlots - 1 % Iteracja po slotach
    % Tworzenie siatki zasobów nośnika dla jednego slotu
    csirsSlotGrid = nrResourceGrid(carrier, csirsPorts); % Generowanie siatki zasobów nośnika

    % Aktualizacja numeru slotu w obiekcie konfiguracji nośnika
    carrier.NSlot = nslot; % 

    csirsInd = nrCSIRSIndices(carrier, csirs); % Generowanie indeksów CSI-RS
    csirsSym = nrCSIRS(carrier, csirs); % Generowanie symboli CSI-RS

    csirsSlotGrid(csirsInd) = csirsSym; % Mapowanie CSI-RS na siatkę slotów

    % Mapowanie portów CSI-RS na anteny nadawcze
    wtx = eye(csirsPorts, nTxAnts); % Generowanie macierzy wag nadawczych  
    txGrid = reshape(reshape(csirsSlotGrid, [], csirsPorts) * wtx, size(csirsSlotGrid, 1), size(csirsSlotGrid, 2), nTxAnts); % Mapowanie portów CSI-RS na anteny nadawcze

    % OFDM modulacja do wygenerowania sygnału w dziedzinie czasu
    txWaveform = nrOFDMModulate(carrier, txGrid); % Modulacja OFDM do wygenerowania sygnału w dziedzinie czasu

    % Dodanie zer na końcu sygnału nadawanego 
    txWaveform = [txWaveform; zeros(maxChDelay, size(txWaveform, 2))]; 

    % Transmisja sygnału przez kanał
    [rxWaveform, pathGains, sampleTimes] = channel(txWaveform); 

    % Generowanie i dodanie szumu AWGN do odbieranego sygnału
    SNR = 10^(SNRdB / 10); % Przeliczenie SNR na wartość liniową
    sigma = 1 / (sqrt(2.0 * nRxAnts * double(waveformInfo.Nfft) * SNR)); % Obliczenie odchylenia standardowego szumu
    noise = sigma * complex(randn(size(rxWaveform)), randn(size(rxWaveform))); % Generowanie szumu
    rxWaveform = rxWaveform + noise; % Dodanie szumu do odbieranego sygnału

    % Praktyczna estymacja przesunięcia czasowego
    [t, mag] = nrTimingEstimate(carrier, rxWaveform, csirsInd, csirsSym); % Estymacja przesunięcia czasowego
    offsetPractical = hSkipWeakTimingOffset(offsetPractical, t, mag); 

    % Pobieranie filtry ścieżek
    pathFilters = getPathFilters(channel); 

    % Perfekcyjna estymacja przesunięcia czasowego
    offsetPerfect = nrPerfectTimingEstimate(pathGains, pathFilters);

    % Korekcja przesunięcia czasowego dla scenariuszy praktycznej i perfekcyjnej estymacji
    rxWaveformPractical = rxWaveform(1 + offsetPractical:end, :); 
    rxWaveformPerfect = rxWaveform(1 + offsetPerfect:end, :); 

    % OFDM demodulacja zsynchronizowanych sygnałów
    rxGridPractical = nrOFDMDemodulate(carrier, rxWaveformPractical); 
    rxGridPerfect = nrOFDMDemodulate(carrier, rxWaveformPerfect);

    % Dodanie zer w przypadku niepełnego slotu po synchronizacji czasowej
    symbPerSlot = carrier.SymbolsPerSlot; % Liczba symboli na slot
    K = size(rxGridPractical, 1); % Rozmiar pierwszego wymiaru siatki odbieranej praktycznie
    LPractical = size(rxGridPractical, 2); % Rozmiar drugiego wymiaru siatki odbieranej praktycznie
    LPerfect = size(rxGridPerfect, 2); % Rozmiar drugiego wymiaru siatki odbieranej perfekcyjnie
    if LPractical < symbPerSlot % Jeśli liczba symboli odbieranych praktycznie jest mniejsza niż wymagana liczba symboli na slot
        rxGridPractical = cat(2, rxGridPractical, zeros(K, symbPerSlot - LPractical, nRxAnts)); % Dodanie zer do siatki praktycznej
    end
    if LPerfect < symbPerSlot % Jeśli liczba symboli odbieranych perfekcyjnie jest mniejsza niż wymagana liczba symboli na slot
        rxGridPerfect = cat(2, rxGridPerfect, zeros(K, symbPerSlot - LPerfect, nRxAnts)); % Dodanie zer do siatki perfekcyjnej
    end
    rxGridPractical = rxGridPractical(:, 1:symbPerSlot, :); % Wybór odpowiedniej liczby symboli dla siatki praktycznej
    rxGridPerfect = rxGridPerfect(:, 1:symbPerSlot, :); % Wybór odpowiedniej liczby symboli dla siatki perfekcyjnej

    % Uwzględnienie tylko symboli i indeksów NZP-CSI-RS do estymacji kanału
    nzpCSIRSSym = csirsSym(csirsSym ~= 0); % Wybór symboli NZP-CSI-RS
    nzpCSIRSInd = csirsInd(csirsSym ~= 0); % Wybór indeksów NZP-CSI-RS

    % Praktyczna estymacja kanału
    [PracticalHest, nVarPractical] = nrChannelEstimate(carrier, rxGridPractical, nzpCSIRSInd, nzpCSIRSSym, 'CDMLengths', cdmLengths, 'AveragingWindow', [0 5]);

    % Perfekcyjna estymacja kanału
    PerfectHest = nrPerfectChannelEstimate(carrier, pathGains, pathFilters, offsetPerfect, sampleTimes); 

    % Perfekcyjna estymacja szumu z realizacji szumu
    noiseGrid = nrOFDMDemodulate(carrier, noise(1 + offsetPerfect:end, :)); % Demodulacja OFDM szumu
    nVarPerfect = var(noiseGrid(:)); % Obliczenie wariancji szumu

    if ~isempty(nzpCSIRSInd) % Jeśli istnieją indeksy NZP-CSI-RS
        % Ustawienie wartości wektora totSlotsBinaryVec dla slotu z NZP-CSI-RS na 1
        totSlotsBinaryVec(nslot + 1) = 1;

        % Estymacja RI na podstawie praktycznej estymacji kanału
        numLayersPractical = hRISelect(carrier, csirs, reportConfig, PracticalHest, nVarPractical, 'MaxSE'); % Wybór liczby warstw na podstawie praktycznej estymacji kanału

        % Estymacja CQI i PMI na podstawie praktycznej estymacji kanału
        [cqiPractical, pmiPractical, cqiInfoPractical, pmiInfoPractical] = hCQISelect(carrier, csirs, reportConfig, numLayersPractical, PracticalHest, nVarPractical); 
        numCodeWordsPr = size(cqiPractical, 2); % Liczba słów kodowych
        numSBs = size(cqiPractical, 1); % Liczba podpasm

        % Przechowywanie wartości CQI, PMI, RI i SINR dla każdego slotu w scenariuszu praktycznej estymacji kanału
        riPracticalPerSlot(1, nslot + 1) = numLayersPractical; % Przechowywanie liczby warstw
        cqiPracticalPerSlot(:, :, nslot + 1) = [cqiPractical NaN(numSBs, 2 - numCodeWordsPr)]; % Przechowywanie CQI
        pmiPracticalPerSlot(nslot + 1) = pmiPractical; % Przechowywanie PMI
        subbandCQIPractical(:, :, nslot + 1) = [cqiInfoPractical.SubbandCQI NaN(numSBs, 2 - numCodeWordsPr)]; % Przechowywanie podpasm CQI
        SINRPerSubbandPerCWPractical(:, :, nslot + 1) = [cqiInfoPractical.SINRPerSubbandPerCW NaN(numSBs, 2 - numCodeWordsPr)]; % Przechowywanie SINR

        % Estymacja RI na podstawie perfekcyjnej estymacji kanału
        numLayersPerfect = hRISelect(carrier, csirs, reportConfig, PerfectHest, nVarPerfect, 'MaxSE'); % Wybór liczby warstw na podstawie perfekcyjnej estymacji kanału

        % Estymacja CQI i PMI na podstawie perfekcyjnej estymacji kanału
        [cqiPerfect, pmiPerfect, cqiInfoPerfect, pmiInfoPerfect] = hCQISelect(carrier, csirs, reportConfig, numLayersPerfect, PerfectHest, nVarPerfect); 
        numCodeWordsPe = size(cqiPerfect, 2); 

        % Przechowywanie wartości CQI, PMI, RI i SINR dla każdego slotu w scenariuszu perfekcyjnej estymacji kanału
        riPerfectPerSlot(1, nslot + 1) = numLayersPerfect; % Przechowywanie liczby warstw
        cqiPerfectPerSlot(:, :, nslot + 1) = [cqiPerfect NaN(numSBs, 2 - numCodeWordsPe)]; % Przechowywanie CQI
        subbandCQIPerfect(:, :, nslot + 1) = [cqiInfoPerfect.SubbandCQI NaN(numSBs, 2 - numCodeWordsPe)]; % Przechowywanie podpasm CQI
        pmiPerfectPerSlot(nslot + 1) = pmiPerfect; % Przechowywanie PMI
        SINRPerSubbandPerCWPerfect(:, :, nslot + 1) = [cqiInfoPerfect.SINRPerSubbandPerCW NaN(numSBs, 2 - numCodeWordsPe)]; % Przechowywanie SINR
    end
end


%% 

% Aktywne numery slotów z obecnym NZP-CSI-RS
activeSlotNum = find(totSlotsBinaryVec); % Wyszukanie aktywnych numerów slotów

% Uzupełnianie zmiennych CQI, PMI i RI wartościami NaN w slotach bez NZP-CSI-RS
[cqiPracticalPerSlot, subbandCQIPractical, pmiPracticalPerSlot, SINRPerSubbandPerCWPractical, cqiPerfectPerSlot, subbandCQIPerfect, pmiPerfectPerSlot, SINRPerSubbandPerCWPerfect, riPracticalPerSlot, riPerfectPerSlot] = fillInactiveSlots(cqiPracticalPerSlot, subbandCQIPractical, pmiPracticalPerSlot, SINRPerSubbandPerCWPractical, cqiPerfectPerSlot, subbandCQIPerfect, pmiPerfectPerSlot, SINRPerSubbandPerCWPerfect, riPracticalPerSlot, riPerfectPerSlot, reportConfig, totSlotsBinaryVec, activeSlotNum); % Uzupełnianie zmiennych wartościami NaN

% Rysowanie wykresów CQI i SINR w szerokopasmie
plotWidebandCQIAndSINR(cqiPracticalPerSlot, cqiPerfectPerSlot, SINRPerSubbandPerCWPractical, SINRPerSubbandPerCWPerfect, activeSlotNum); % Rysowanie wykresów CQI i SINR w szerokopasmie

% Numer slotu do rysowania wykresów CQI i SINR w podpasie
slotNumForSBType1 = 0; % Ustalenie numeru slotu na 0 (liczone od zera)
plotSubbandCQIAndSINR(subbandCQIPractical, subbandCQIPerfect, SINRPerSubbandPerCWPractical, SINRPerSubbandPerCWPerfect, activeSlotNum, slotNumForSBType1); % Rysowanie wykresów CQI i SINR w podpasie

% Flaga do wyświetlania wykresów PMI i RI dla codebooku typu I
showType1PMIandRIPlots = true; % Ustawienie flagi do wyświetlania wykresów PMI i RI dla codebooku typu I

% Flaga do wyświetlania wykresów PMI i RI dla codebooku typu II
showType2PMIRI = true; % Ustawienie flagi do wyświetlania wykresów PMI i RI dla codebooku typu II

% Rysowanie wykresów PMI i RI dla codebooku typu I
if (strcmpi(reportConfig.CodebookType, 'Type1SinglePanel') || strcmpi(reportConfig.CodebookType, 'Type1MultiPanel')) && showType1PMIandRIPlots % Sprawdzenie warunków do wyświetlania wykresów PMI i RI dla codebooku typu I
    plotType1PMIAndRI(pmiPracticalPerSlot, pmiPerfectPerSlot, riPracticalPerSlot, riPerfectPerSlot, activeSlotNum, slotNumForSBType1); % Rysowanie wykresów PMI i RI dla codebooku typu I

% Rysowanie wykresów PMI i RI dla codebooku typu II
    slotNumForType2PMI = 0; % Ustawienie numeru slotu na 0 (liczone od zera)
    plotType2PMIAndRI(pmiPracticalPerSlot, pmiPerfectPerSlot, riPracticalPerSlot, riPerfectPerSlot, reportConfig.PanelDimensions, reportConfig.NumberOfBeams, activeSlotNum, slotNumForType2PMI); % Rysowanie wykresów PMI i RI dla codebooku typu II
end







%%%############################################# FUNKCJE ###############%%

% definicje funkcji 
% Sprawdza konfigurację CSI-RS dla danego nośnika.
% %validateCSIRSConfig sprawdza liczbę portów, długości CDM, liczbę anten i indeksy CSI-RS.

% Zwraca długości CDM na podstawie konfiguracji CSI-RS.
% %getCDMLengths pobiera typ CDM i zwraca odpowiednie długości.

% Uzupełnia zmienne CQI, PMI i RI wartościami NaN w nieobecnych slotach.
% %fillInactiveSlots oblicza indeksy brakujących slotów i uzupełnia zmienne.

% Tworzy wykresy CQI i SINR w paśmie dla każdego kodu w aktywnych slotach.
% %plotWidebandCQIAndSINR sprawdza dane CQI i tworzy wykresy.

% Tworzy wykresy SINR i CQI dla każdego kodu na podpasmach w określonym slocie.
% %plotSubbandCQIAndSINR sprawdza dane CQI i tworzy wykresy na podpasmach.

% Tworzy wykresy PMI i RI dla każdego podpasma lub kodu w określonym slocie.
% %plotType1PMIAndRI sprawdza dane PMI i RI i tworzy wykresy.

% Tworzy siatkę wiązek/promieni/beams, podkreślając używane do generowania macierzy precodowania.
% %plotType2PMIAndRI wyszukuje i rysuje odpowiednie wiązek/promieni/beams na siatce.


function validateCSIRSConfig(carrier,csirs,nTxAnts)
%   Validates the CSI-RS configuration, given the carrier specific
%   configuration object, CSI-RS configuration object, and the number of
%   transmit antennas.

    % Validate the number of CSI-RS ports
    if ~isscalar(unique(csirs.NumCSIRSPorts))
        error('nr5g:InvalidCSIRSPorts', ...
            'All the CSI-RS resources must be configured to have the same number of CSI-RS ports.');
    end

    % Validate the CDM lengths
    if ~iscell(csirs.CDMType)
        cdmType = {csirs.CDMType};
    else
        cdmType = csirs.CDMType;
    end
    if (~all(strcmpi(cdmType,cdmType{1})))
        error('nr5g:InvalidCSIRSCDMTypes', ...
            'All the CSI-RS resources must be configured to have the same CDM lengths.');
    end
    if nTxAnts ~= csirs.NumCSIRSPorts(1)
        error('nr5g:InvalidNumTxAnts',['Number of transmit antennas (' num2str(nTxAnts) ...
            ') must be equal to the number of CSI-RS ports (' num2str(csirs.NumCSIRSPorts(1)) ').']);
    end

    % Check for the overlap between the CSI-RS indices
    csirsInd = nrCSIRSIndices(carrier,csirs,"OutputResourceFormat",'cell');
    numRes = numel(csirsInd);
    csirsIndAll = cell(1,numRes);
    ratioVal = csirs.NumCSIRSPorts(1)/prod(getCDMLengths(csirs));
    for resIdx = 1:numRes
        if ~isempty(csirsInd{resIdx})
            grid = nrResourceGrid(carrier,csirs.NumCSIRSPorts(1));
            [~,tempInd] = nrExtractResources(csirsInd{resIdx},grid);
            if numel(tempInd)/numel(csirsInd{resIdx}) ~= ratioVal
                error('nr5g:OverlappedCSIRSREsSingleResource',['CSI-RS indices of resource ' ...
                    num2str(resIdx) ' must be unique. Try changing the symbol or subcarrier locations.']);
            end
            csirsIndAll{resIdx} = tempInd(:);
            for idx = 1:resIdx-1
                overlappedInd = ismember(csirsIndAll{idx},csirsIndAll{resIdx});
                if any(overlappedInd)
                    error('nr5g:OverlappedCSIRSREsMultipleResources',['The resource elements of the ' ...
                        'configured CSI-RS resources must not overlap. Try changing the symbol or ' ...
                        'subcarrier locations of CSI-RS resource ' num2str(idx) ' and resource ' num2str(resIdx) '.']);
                end
            end
        end
    end
end

function cdmLengths = getCDMLengths(csirs)
%   Returns the CDM lengths, given the CSI-RS configuration object.

    CDMType = csirs.CDMType;
    if ~iscell(csirs.CDMType)
        CDMType = {csirs.CDMType};
    end
    CDMTypeOpts = {'noCDM','fd-CDM2','CDM4','CDM8'};
    CDMLengthOpts = {[1 1],[2 1],[2 2],[2 4]};
    cdmLengths = CDMLengthOpts{strcmpi(CDMTypeOpts,CDMType{1})};
end

function [cqiPracticalPerSlot,subbandCQIPractical,pmiPracticalPerSlot,SINRPerSubbandPerCWPractical,cqiPerfectPerSlot, ...
    subbandCQIPerfect,pmiPerfectPerSlot,SINRPerSubbandPerCWPerfect,riPracticalPerSlot,riPerfectPerSlot] = fillInactiveSlots(cqiPracticalPerSlot, ...
    subbandCQIPractical,pmiPracticalPerSlot,SINRPerSubbandPerCWPractical,cqiPerfectPerSlot,subbandCQIPerfect,pmiPerfectPerSlot, ...
    SINRPerSubbandPerCWPerfect,riPracticalPerSlot,riPerfectPerSlot,reportConfig,totSlotsBinaryVec,activeSlots)
%   Returns the CQI, PMI, and RI related variables filled with NaNs in the
%   slots where NZP-CSI-RS is not present according to the codebook type from
%   the report configuration structure. Note that the CQI, PMI, and RI
%   variables are returned as empty if there are no NZP-CSI-RS resources,
%   that is, no active slots in the entire simulation duration.

    % Compute the indices of the slots and the number of slots in which
    % NZP-CSI-RS is not present
    inactiveSlotIdx = ~totSlotsBinaryVec;
    numInactiveSlots = nnz(inactiveSlotIdx);
    
    if ~isempty(activeSlots)
        numCQISBs = size(cqiPracticalPerSlot,1);
    
        % Get the codebook type
        codebookType = 'Type1SinglePanel';
        if isfield(reportConfig,'CodebookType')
            codebookType = validatestring(reportConfig.CodebookType,{'Type1SinglePanel','Type1MultiPanel','Type2','eType2'},'fillInactiveSlots','CodebookType field');
        end
    
        % Fill the CQI, PMI, and RI variables with NaNs in the slots where NZP-CSI-RS is
        % not present
        cqiPracticalPerSlot(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        subbandCQIPractical(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        SINRPerSubbandPerCWPractical(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        cqiPerfectPerSlot(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        subbandCQIPerfect(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        SINRPerSubbandPerCWPerfect(:,:,inactiveSlotIdx) = NaN(numCQISBs,2,numInactiveSlots);
        riPracticalPerSlot(inactiveSlotIdx) = NaN;
        riPerfectPerSlot(inactiveSlotIdx) = NaN;
    
        numi1Indices = 3;
        numi2Indices = 1;
        if strcmpi(codebookType,'Type1MultiPanel')
            numi1Indices = 6;
            numi2Indices = 3;
        end
        numPMISBs = size(pmiPerfectPerSlot(activeSlots(1)).i2,2);
        [pmiPerfectPerSlot(inactiveSlotIdx),pmiPracticalPerSlot(inactiveSlotIdx)] = deal(struct('i1',NaN(1,numi1Indices),'i2',NaN(numi2Indices,numPMISBs)));
    end
end

function plotWidebandCQIAndSINR(cqiPracticalPerSlot,cqiPerfectPerSlot,SINRPerSubbandPerCWPractical,SINRPerSubbandPerCWPerfect,activeSlotNum)
%   Plots the wideband SINR and wideband CQI values for each codeword
%   across all specified active slots (1-based) (in which the CQI is
%   reported as other than NaN) for practical and perfect channel
%   estimation cases.

    % Check if there are no slots in which NZP-CSI-RS is present
    if isempty(activeSlotNum)
        disp('No CQI data to plot, because there are no slots in which NZP-CSI-RS is present.');
        return;
    end
    cqiPracticalPerCW = permute(cqiPracticalPerSlot(1,:,:),[1 3 2]);
    cqiPerfectPerCW = permute(cqiPerfectPerSlot(1,:,:),[1 3 2]);
    SINRPerCWPractical = permute(SINRPerSubbandPerCWPractical(1,:,:),[1 3 2]);
    SINRPerCWPerfect = permute(SINRPerSubbandPerCWPerfect(1,:,:),[1 3 2]);

    % Extract wideband CQI indices for slots where NZP-CSI-RS is present
    cqiPracticalPerCWActiveSlots = cqiPracticalPerCW(1,activeSlotNum,:);
    cqiPerfectPerCWActiveSlots = cqiPerfectPerCW(1,activeSlotNum,:);
    widebandSINRPractical = 10*log10(SINRPerCWPractical(1,activeSlotNum,:));
    widebandSINRPerfect = 10*log10(SINRPerCWPerfect(1,activeSlotNum,:));

    if isempty(reshape(cqiPracticalPerCWActiveSlots(:,:,1),1,[]))
        disp('No CQI data to plot, because all CQI values are NaNs.');
        return;
    end

    figure();
    plotWBCQISINR(widebandSINRPerfect,widebandSINRPractical,211,activeSlotNum,'SINR');
    plotWBCQISINR(cqiPerfectPerCWActiveSlots,cqiPracticalPerCWActiveSlots,212,activeSlotNum,'CQI');
end

function plotWBCQISINR(perfectVals,practicalVals,subplotIdx,activeSlotNum,inpText)
%   Plots the wideband SINR and wideband CQI values for each codeword
%   across all specified active slots (1-based) (in which the CQI is
%   reported as other than NaN) for practical and perfect channel
%   estimation cases.

    subplot(subplotIdx)
    plot(perfectVals(:,:,1),'r-o');
    hold on;
    plot(practicalVals(:,:,1),'b-*');
    if ~all(isnan(perfectVals(:,:,2))) % Two codewords
        hold on;
        plot(perfectVals(:,:,2),'r:s');
        hold on;
        plot(practicalVals(:,:,2),'b:d');
        title(['Wideband ' inpText ' Values for Codeword 1&2']);
        legend({'Codeword 1:Perfect channel est.','Codeword 1:Practical channel est.','Codeword 2:Perfect channel est.','Codeword 2:Practical channel est.'});
    else
        title(['Wideband ' inpText ' Values for Codeword 1']);
        legend({'Codeword 1:Perfect channel est.','Codeword 1:Practical channel est.'});
    end
    xlabel('Slots');
    if strcmpi(inpText,'SINR')
        units = ' in dB';
    else
        units = '';
    end
    ylabel(['Wideband ' inpText ' Values' units]);
    xticks(1:size(perfectVals,2));
    xTickLables = num2cell(activeSlotNum(:)-1);
    xticklabels(xTickLables);
    [lowerBound,upperBound] = bounds([practicalVals(:);perfectVals(:)]);
    ylim([lowerBound-1 upperBound+3.5]);
end

function plotSubbandCQIAndSINR(subbandCQIPractical,subbandCQIPerfect,SINRPerCWPractical,SINRPerCWPerfect,activeSlotNum,nslot)
%   Plots the SINR and CQI values for each codeword across all the subbands
%   for practical and perfect channel estimation cases for the given slot
%   number (0-based) among all specified active slots (1-based). The
%   function does not plot the values if CQIMode is 'Wideband' or if the
%   CQI and SINR values are all NaNs in the given slot.

    % Check if there are no slots in which NZP-CSI-RS is present
    if isempty(activeSlotNum)
        disp('No CQI data to plot, because there are no slots in which NZP-CSI-RS is present.');
        return;
    end
    numSubbands = size(subbandCQIPractical,1);
    if numSubbands > 1 && ~any(nslot+1 == activeSlotNum) % Check if the CQI values are reported in the specified slot
        disp(['For the specified slot (' num2str(nslot) '), CQI values are not reported. Please choose another slot number.']);
        return;
    end

    % Plot subband CQI values
    if numSubbands > 1 % Subband mode
        subbandCQIPerCWPractical = subbandCQIPractical(2:end,:,nslot+1);
        subbandCQIPerCWPerfect = subbandCQIPerfect(2:end,:,nslot+1);
        subbandSINRPerCWPractical = 10*log10(SINRPerCWPractical(2:end,:,nslot+1));
        subbandSINRPerCWPerfect = 10*log10(SINRPerCWPerfect(2:end,:,nslot+1));
        figure();
        plotSBCQISINR(subbandSINRPerCWPerfect,subbandSINRPerCWPractical,numSubbands,211,nslot,'SINR')
        plotSBCQISINR(subbandCQIPerCWPerfect,subbandCQIPerCWPractical,numSubbands,212,nslot,'CQI');
    end
end

function plotSBCQISINR(perfectVals,practicalVals,numSubbands,subplotIdx,nslot,inpText)
%   Plots the SINR and CQI values for each codeword across all the subbands
%   for practical and perfect channel estimation cases for the given slot
%   number (0-based). The function does not plot the values if CQIMode is
%   'Wideband' or if the CQI and SINR values are all NaNs in the given
%   slot.

    subplot(subplotIdx)
    plot(perfectVals(:,1),'ro-');
    hold on;
    plot(practicalVals(:,1),'b*-');
    if ~all(isnan(perfectVals(:,2))) % Two codewords
        hold on;
        plot(perfectVals(:,2),'rs:');
        hold on;
        plot(practicalVals(:,2),'bd:');
        legend({'Codeword 1:Perfect channel est.','Codeword 1:Practical channel est.','Codeword 2:Perfect channel est.','Codeword 2:Practical channel est.'});
        title(['Estimated Subband ' inpText ' Values for Codeword 1&2 in Slot ' num2str(nslot)]);
    else % Single codeword
        legend({'Codeword 1:Perfect channel est.','Codeword 1:Practical channel est.'});
        title(['Estimated Subband ' inpText ' Values for Codeword 1 in Slot ' num2str(nslot)]);
    end

    if strcmpi(inpText,'SINR')
        units = ' in dB';
    else
        units = '';
    end
    xlabel('Subbands');
    ylabel(['Subband ' inpText ' Values' units]);
    xticks(1:numSubbands);
    xTickLables = num2cell(1:numSubbands);
    xticklabels(xTickLables);
    xlim([0 numSubbands+1]);
    [lowerBound,upperBound] = bounds([perfectVals(:);practicalVals(:)]);
    ylim([lowerBound-1 upperBound+3.5]);
end

function plotType1PMIAndRI(pmiPracticalPerSlot,pmiPerfectPerSlot,riPracticalPerSlot,riPerfectPerSlot,activeSlotNum,nslot)
%   Plots the RI and PMI i1 indices across all specified active slots
%   (1-based), for practical and perfect channel estimation scenarios. The
%   function also plots the i2 indices of practical and perfect channel
%   estimation scenarios across all specified active slots when the PMI
%   mode is 'Wideband' or plots i2 indices across all the subbands for the
%   specified slot number (0-based) when the PMI mode is 'Subband'.

    % Check if there are no slots in which NZP-CSI-RS is present
    if isempty(activeSlotNum)
        disp('No PMI and RI data to plot, because there are no slots in which NZP-CSI-RS is present.');
        return;
    end
    
    numi1Indices = numel(pmiPracticalPerSlot(activeSlotNum(1)).i1);
    if numi1Indices == 6
        codebookType = 'Type1MultiPanel';
    else
        codebookType = 'Type1SinglePanel';
    end
    
    % Extract wideband PMI indices (i1 values) for slots where NZP-CSI-RS
    % is present
    i1PerfectValsActiveSlots = reshape([pmiPerfectPerSlot(activeSlotNum).i1],numi1Indices,[])';
    i1PracticalValsActiveSlots = reshape([pmiPracticalPerSlot(activeSlotNum).i1],numi1Indices,[])';
    
    if isempty(i1PerfectValsActiveSlots)
        disp('No PMI and RI data to plot, because all PMI and RI values are NaNs.');
        return;
    end
    
    figure;
    % Plot RI
    plotRI(riPracticalPerSlot,riPerfectPerSlot,activeSlotNum,411);
    
    % Extract and plot i11 indices
    i11PerfectVals = i1PerfectValsActiveSlots(:,1);
    i11PracticalVals = i1PracticalValsActiveSlots(:,1);
    plotIxxIndices(i11PerfectVals,i11PracticalVals,activeSlotNum,412,'i11');

    % Extract and plot i12 indices
    i12PerfectVals = i1PerfectValsActiveSlots(:,2);
    i12PracticalVals = i1PracticalValsActiveSlots(:,2);
    plotIxxIndices(i12PerfectVals,i12PracticalVals,activeSlotNum,413,'i12');

    % Extract and plot i13 indices
    i13PerfectVals = i1PerfectValsActiveSlots(:,3);
    i13PracticalVals = i1PracticalValsActiveSlots(:,3);
    plotIxxIndices(i13PerfectVals,i13PracticalVals,activeSlotNum,414,'i13');
    
    % Plot the i141, i142 and i143 indices in type I multi-panel case
    if strcmpi(codebookType,'Type1MultiPanel')
        figure()
        % Extract and plot i141 indices
        i141PerfectVals = i1PerfectValsActiveSlots(:,4);
        i141PracticalVals = i1PracticalValsActiveSlots(:,4);
        plotIxxIndices(i141PerfectVals,i141PracticalVals,activeSlotNum,311,'i141');

        % Extract and plot i142 indices
        i142PerfectVals = i1PerfectValsActiveSlots(:,5);
        i142PracticalVals = i1PracticalValsActiveSlots(:,5);
        plotIxxIndices(i142PerfectVals,i142PracticalVals,activeSlotNum,312,'i142');
    
        % Extract and plot i143 indices
        i143PerfectVals = i1PerfectValsActiveSlots(:,6);
        i143PracticalVals = i1PracticalValsActiveSlots(:,6);
        plotIxxIndices(i143PerfectVals,i143PracticalVals,activeSlotNum,313,'i143');
    end

    % Get the number of subbands
    numSubbands = size(pmiPracticalPerSlot(activeSlotNum(1)).i2,2);
    % Get the number of i2 indices according to codebook type
    numi2Indices = 1;
    if strcmpi(codebookType,'Type1MultiPanel')
        numi2Indices = 3;
    end

    % Get number of active slots
    numActiveSlots = numel(activeSlotNum);
    % Extract i2 values
    i2PerfectVals = reshape([pmiPerfectPerSlot(activeSlotNum).i2],[numSubbands,numi2Indices,numActiveSlots]);     % Of size numActiveSlots-by-numi2Indices-numSubbands
    i2PracticalVals = reshape([pmiPracticalPerSlot(activeSlotNum).i2],[numSubbands,numi2Indices,numActiveSlots]); % Of size numActiveSlots-by-numi2Indices-numSubbands

    % Plot i2 values
    if numSubbands == 1 % Wideband mode
        figure;

        % In type I single-panel case, there is only one i2 index. The
        % first column of i2PerfectVals and i2PracticalVals corresponds to
        % i2 index. In type I multi-panel case, the i2 values are a set of
        % three indices i20, i21, and i22. Each column of i2PerfectVals and
        % i2PracticalVals correspond to i20, i21, and i22 indices. Extract
        % and plot the respective index values
        if strcmpi(codebookType,'Type1SinglePanel')
            % Extract and plot i2 values in each slot
            i2PerfectVals = reshape(i2PerfectVals(:,1,:),[],numActiveSlots).';
            i2PracticalVals = reshape(i2PracticalVals(:,1,:),[],numActiveSlots).';
            plotIxxIndices(i2PerfectVals,i2PracticalVals,activeSlotNum,111,'i2');
        else
            % Extract and plot i20 values in each slot
            i20PerfectVals = reshape(i2PerfectVals(:,1,:),[],numActiveSlots).';
            i20PracticalVals = reshape(i2PracticalVals(:,1,:),[],numActiveSlots).';
            plotIxxIndices(i20PerfectVals,i20PracticalVals,activeSlotNum,311,'i20');

            % Extract and plot i21 values in each slot
            i21PerfectVals = reshape(i2PerfectVals(:,2,:),[],numActiveSlots).';
            i21PracticalVals = reshape(i2PracticalVals(:,2,:),[],numActiveSlots).';
            plotIxxIndices(i21PerfectVals,i21PracticalVals,activeSlotNum,312,'i21');

            % Extract and plot i22 values in each slot
            i22PerfectVals = reshape(i2PerfectVals(:,3,:),[],numActiveSlots).';
            i22PracticalVals = reshape(i2PracticalVals(:,3,:),[],numActiveSlots).';
            plotIxxIndices(i22PerfectVals,i22PracticalVals,activeSlotNum,313,'i22');
        end
    else % Subband mode
        if any(nslot+1 == activeSlotNum)
    
            % In subband mode, plot the PMI i2 indices corresponding to the
            % specified slot number
            figure;

            if strcmpi(codebookType,'Type1SinglePanel')
                % Extract and plot i2 values
                pmiSBi2Perfect = pmiPerfectPerSlot(nslot+1).i2(1,:);
                pmiSBi2Practical = pmiPracticalPerSlot(nslot+1).i2(1,:);
                plotI2xIndices_SB(pmiSBi2Perfect,pmiSBi2Practical,numSubbands,nslot,111,'i2');
            else
                % Extract and plot i20 values
                pmiSBi20Perfect = pmiPerfectPerSlot(nslot+1).i2(1,:);
                pmiSBi20Practical = pmiPracticalPerSlot(nslot+1).i2(1,:);
                plotI2xIndices_SB(pmiSBi20Perfect,pmiSBi20Practical,numSubbands,nslot,311,'i20');
                
                % Extract and plot i21 values
                pmiSBi21Perfect = pmiPerfectPerSlot(nslot+1).i2(2,:);
                pmiSBi21Practical = pmiPracticalPerSlot(nslot+1).i2(2,:);
                plotI2xIndices_SB(pmiSBi21Perfect,pmiSBi21Practical,numSubbands,nslot,312,'i21');
    
                % Extract and plot i22 values
                pmiSBi22Perfect = pmiPerfectPerSlot(nslot+1).i2(3,:);
                pmiSBi22Practical = pmiPracticalPerSlot(nslot+1).i2(3,:);
                plotI2xIndices_SB(pmiSBi22Perfect,pmiSBi22Practical,numSubbands,nslot,313,'i22');
            end
        else
            disp(['For the specified slot (' num2str(nslot) '), PMI i2 indices are not reported. Please choose another slot number.'])
        end
    end
end

function plotType2PMIAndRI(pmiPracticalPerSlot,pmiPerfectPerSlot,riPracticalPerSlot,riPerfectPerSlot,panelDims,numBeams,activeSlotNum,nslot)
%   Plots the grid of beams by highlighting the beams that are used for the
%   precoding matrix generation for the specified slot number (0-based),
%   for practical and perfect channel estimation scenarios.

    % Check if there are no slots in which NZP-CSI-RS is present
    if isempty(activeSlotNum)
        disp('No PMI and RI data to plot, because there are no slots in which NZP-CSI-RS is present.');
        return;
    end
    plotRI(riPracticalPerSlot,riPerfectPerSlot,activeSlotNum,111);
    if ~any(nslot+1 == activeSlotNum)
        disp(['For the specified slot (' num2str(nslot) '), PMI values are not reported. Please choose another slot number.']);
    else
        pmiPractical = pmiPracticalPerSlot(nslot+1);
        pmiPerfect = pmiPerfectPerSlot(nslot+1);
        figure();
        plotType2GridOfBeams(pmiPractical,panelDims,numBeams,'Practical Channel Estimation Scenario',1);
        hold on;
        plotType2GridOfBeams(pmiPerfect,panelDims,numBeams,'Perfect Channel Estimation Scenario',2);
    end
end

function plotRI(riPracticalPerSlot,riPerfectPerSlot,activeSlotNum,subplotIndex)
%   Plots the RI values across all specified active slots (1-based), for
%   practical and perfect channel estimation scenarios.

    % Get number of active slots
    numActiveSlots = numel(activeSlotNum);

    % Extract RI values for slots where NZP-CSI-RS is present
    RIPerfectValsActiveSlots = riPerfectPerSlot(activeSlotNum)';
    RIPracticalValsActiveSlots = riPracticalPerSlot(activeSlotNum)';
    
    if isempty(RIPerfectValsActiveSlots)
        disp('No RI data to plot, because all RI values are NaNs.');
        return;
    end
    
    figure;
    subplot(subplotIndex);
    plot(RIPerfectValsActiveSlots,'r-o');
    hold on;
    plot(RIPracticalValsActiveSlots,'b-*');
    xlabel('Slots')
    ylabel('RI Values');
    xticks(1:numActiveSlots);
    xTickLables = num2cell(activeSlotNum(:)-1);
    xticklabels(xTickLables);
    [~,upperBound] = bounds([RIPerfectValsActiveSlots; RIPracticalValsActiveSlots]);
    xlim([0 numActiveSlots+8]);
    ylim([0 upperBound+1]);
    yticks(0:upperBound+1);
    title('RI Values')
    legend({'Perfect channel est.','Practical channel est.'});
end

function plotType2GridOfBeams(PMISet,panelDims,numBeams,chEstType,subplotNum)
%   Plots the grid of beams by highlighting the beams that are used for the
%   type II codebook based precoding matrix generation.

    N1 = panelDims(1);
    N2 = panelDims(2);    
    % Get the oversampling factors
    O1 = 4;
    O2 = 1 + 3*(N2 ~= 1);

    % Extract q1, q2 values
    qSet = PMISet.i1(1:2);
    q1 = qSet(1)-1;
    q2 = qSet(2)-1;

    % Extract i12 value
    i12 = PMISet.i1(3);
    s = 0;
    % Find the n1, n2 values for all the beams, as defined in TS 38.214
    % Section 5.2.2.2.3
    n1_i12 = zeros(1,numBeams);
    n2_i12 = zeros(1,numBeams);
    for beamIdxI = 0:numBeams-1
        i12minussVal = i12 - s;
        xValues = numBeams-1-beamIdxI:N1*N2-1-beamIdxI;
        CValues = zeros(numel(xValues),1);
        for xIdx = 1:numel(xValues)
            if xValues(xIdx) >= numBeams-beamIdxI
                CValues(xIdx) = nchoosek(xValues(xIdx),numBeams-beamIdxI);
            end
        end
        indices = i12minussVal >= CValues;
        maxIdx = find(indices,1,'last');
        xValue = xValues(maxIdx);
        ei = CValues(maxIdx);
        s = s+ei;
        ni = N1*N2 - 1 - xValue;
        n1_i12(beamIdxI+1) = mod(ni,N1);
        n2_i12(beamIdxI+1) = (ni-n1_i12(beamIdxI+1))/N1;
    end
    m1 = O1*(0:N1-1) + q1;
    m2 = O2*(0:N2-1) + q2;

    % Calculate the indices of orthogonal basis set which corresponds to
    % the reported i12 value
    m1_LBeams = O1*(n1_i12) + q1;
    m2_LBeams = O2*(n2_i12) + q2;
    OrthogonalBeams = [repmat(m1,1,length(m2));reshape(repmat(m2,length(m1),1),1,[])]';

    % Plot the grid of beams
    numCirlcesInRow = N1*O1;
    numCirlcesInCol = N2*O2;
    subplot(2,1,subplotNum);
    circleRadius = 1;
    for colIdx = 0:numCirlcesInCol-1
        for rowIdx = 0:numCirlcesInRow-1
            p = nsidedpoly(1000, 'Center', [2*rowIdx 2*colIdx], 'Radius', circleRadius);
            if any(prod(OrthogonalBeams == [rowIdx colIdx],2))
                h2 = plot(p, 'FaceColor', 'w','EdgeColor','r','LineWidth',2.5);
                hold on;
                if any(prod([m1_LBeams' m2_LBeams'] == [rowIdx colIdx],2))
                    h3 = plot(p, 'FaceColor', 'g','LineStyle','-.');                
                end
            else
                h1 = plot(p, 'FaceColor', 'w');
            end
            hold on;
        end
    end
    rowLength = 2*circleRadius*O1;
    colLength = 2*circleRadius*O2;
    for n2 = 0:N2-1
        for n1 = 0:N1-1
            x1 = -1*circleRadius + rowLength*n1;
            x2 = x1 + rowLength;
            y1 = -1*circleRadius + colLength*n2;
            y2 = y1 + colLength;
            x = [x1, x2, x2, x1, x1];
            y = [y1, y1, y2, y2, y1];
            plot(x, y, 'b-', 'LineWidth', 2);
            hold on;
        end
    end
    
    xlabel('N1O1 beams');
    ylabel('N2O2 beams');
    axis equal;
    set(gca,'xtick',[],'ytick',[]);
    legend([h1 h2 h3],{'Oversampled DFT beams',['Orthogonal basis set with [q1 q2] = [' num2str(q1) ' ' num2str(q2) ']'],'Selected beam group'},'Location','northeast');
    title(['Grid of Beams or DFT Vectors for ' chEstType]);
end

function plotIxxIndices(ixxPerfectVals,ixxPracticalVals,activeSlotNum,subplotInp,pmiIdxType)
%   Plots i11, i12, i13 indices in case of type I single-panel codebooks
%   and plots i141, i142, and i143 in case of type I multi-panel codebooks.

    % Plot ixx values
    subplot(subplotInp)
    plot(ixxPerfectVals,'r-o');
    hold on;
    plot(ixxPracticalVals,'b-*');
    xlabel('Slots')
    ylabel([pmiIdxType ' Indices']);
    % Get number of active slots
    numActiveSlots = numel(activeSlotNum);
    xticks(1:numActiveSlots);
    xTickLables = num2cell(activeSlotNum(:)-1);
    xticklabels(xTickLables);
    [lowerBound,upperBound] = bounds([ixxPerfectVals; ixxPracticalVals]);
    xlim([0 numActiveSlots+8]);
    ylim([lowerBound-2 upperBound+2]);
    title(['PMI: ' pmiIdxType ' Indices']);
    legend({'Perfect channel est.','Practical channel est.'});
end

function plotI2xIndices_SB(pmiSBi2Perfect,pmiSBi2Practical,numSubbands,nslot,subplotInp,pmiIdxType)
%   Plots i2 indices in case of type I single-panel codebooks and plots
%   i20, i21, and i22 in case of type I multi-panel codebooks.

    subplot(subplotInp)
    plot(pmiSBi2Perfect,'r-o');
    hold on;
    plot(pmiSBi2Practical,'b-*');
    title(['PMI: ' pmiIdxType ' Indices for All Subbands in Slot ' num2str(nslot)]);
    xlabel('Subbands')
    ylabel([pmiIdxType ' Indices']);
    xticks(1:numSubbands);
    xticklabels(num2cell(1:numSubbands));
    [lowerBound,upperBound] = bounds([pmiSBi2Perfect pmiSBi2Practical]);
    yticks(lowerBound:upperBound);
    yticklabels(num2cell(lowerBound:upperBound));
    xlim([0 numSubbands+1])
    ylim([lowerBound-1 upperBound+1]);
    legend({'Perfect channel est.','Practical channel est.'});
end