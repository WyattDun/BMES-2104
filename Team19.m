% Originally created by Luke Riexinger March 2020
% Edited and modified by Sara Arena April 1 2022
% Edited and modified by Wyatt Dunbar on November 12 2024
clear
close all
clc
%% Simulation Time and Size
areaSize =  150;
total_time =  1000; %frames
dt =  1; %time step
%% Healthy Cell Properties
CellNum =  25; %number of healthy cells at start of simulation
ContactThreshold = 4; % number of contacts required for a healthy cell to be removed
CellLocation = zeros([CellNum, 2]); % Initialize location of healthy cells with zeros
ContactCount = zeros([CellNum, 1]); % Initialize contact count for each healthy cell
for cell = 1 : CellNum
   CellLocation(cell, :) = PlaceCell(areaSize); % place healthy cells in random location within the area
end
%% Immunocompromised Cell Properties
ImmunoCellNum =  10; %number of immunocompromised cells at start of simulation
ImmunoCellLocation = zeros([ImmunoCellNum, 2]); % Initialize location of immunocompromised cells with zeros
for cell = 1 : ImmunoCellNum
   ImmunoCellLocation(cell, :) = PlaceCell(areaSize); % place immunocompromised cells in random location within the area
end
%% Virus Properties
virusNum =  10; %number of viruses at start of simulation
virusLocation = zeros([virusNum, 2]);
for virus = 1 : virusNum
   virusLocation(virus, :) = PlaceCell(areaSize); % place viruses in random location within the area
end
virusTarget = randi([1, CellNum + ImmunoCellNum], 1, virusNum); % randomly select which healthy or immunocompromised cell each virus will attack
virusesGenerated = 1; %the number of new viruses that are created from an infected cell
virusSpeed = 1; %how fast the virus can move towards the targeted cell
%% Immune Cell Properties
ImmuneNum =  2; % number of immune cells at start of simulation
immuneLocation = zeros([ImmuneNum, 2]); % Initialize location of immune cells with zeros
for immune = 1 : ImmuneNum
   immuneLocation(immune, :) = PlaceCell(areaSize); %randomly place immune cells within area
end
immuneSpeed = 1.3; %how fast immune cells can move
immuneReplicateNum = 20; %how many viruses must be eaten before immune cell calls for help
immuneVirusesConsumed = ones([ImmuneNum, 1]); %variable we will use to keep track of how many viruses each immune cell eats
%% Plot the initial locations of all agents
figure()
plot(CellLocation(:, 1), CellLocation(:, 2), 'ok')
hold on
plot(ImmunoCellLocation(:, 1), ImmunoCellLocation(:, 2), 'ob')
plot(virusLocation(:, 1), virusLocation(:,2), '.r')
plot(immuneLocation(:, 1), immuneLocation(:,2), 'og')
xlim([-areaSize, areaSize])
ylim([-areaSize, areaSize])
hold off
%% Simulation
for t = 0:dt:total_time % loop through each instant in time based on time step and simulation time
   % Update each virus in simulation
   for virus = 1:virusNum
       % Determine if the target is a healthy cell or an immunocompromised cell
       if virusTarget(virus) <= CellNum
           targetLocation = CellLocation(virusTarget(virus), :);
           targetType = 'Healthy';
       else
           targetLocation = ImmunoCellLocation(virusTarget(virus) - CellNum, :);
           targetType = 'Immunocompromised';
       end
      
       % Move the viruses
       distanceToCell = ComputeDistance(virusLocation(virus, :), targetLocation);
       angleToCell = ComputeAngle(virusLocation(virus , :), targetLocation);
       virusLocation(virus, :) = virusLocation(virus, :) + [cos(angleToCell), sin(angleToCell)] * virusSpeed;
      
       % If the virus reaches the cell, count the contact and handle infection
       if distanceToCell <= 1
           if strcmp(targetType, 'Healthy')
               ContactCount(virusTarget(virus)) = ContactCount(virusTarget(virus)) + 1;
               if ContactCount(virusTarget(virus)) >= ContactThreshold
                   virusNum = virusNum + virusesGenerated;
                   for newVirus = 1: virusesGenerated
                       virusLocation = [virusLocation; PlaceCell(areaSize)]; % Place new viruses at random location
                       virusTarget = [virusTarget, randi([1, CellNum + ImmunoCellNum])]; % Assign a random cell for the new viruses to target
                   end
                   CellLocation(virusTarget(virus), :) = PlaceCell(areaSize); % Remove "dead" healthy cell and place a new healthy cell randomly to replace the dead one
                   ContactCount(virusTarget(virus)) = 0; % Reset the contact count for the new cell
               end
           else % Immunocompromised cell
               virusNum = virusNum + virusesGenerated;
               for newVirus = 1: virusesGenerated
                   virusLocation = [virusLocation; PlaceCell(areaSize)]; % Place new viruses at random location
                   virusTarget = [virusTarget, randi([1, CellNum + ImmunoCellNum])]; % Assign a random cell for the new viruses to target
               end
               ImmunoCellLocation(virusTarget(virus) - CellNum, :) = PlaceCell(areaSize); % Remove "dead" immunocompromised cell and place a new one randomly
           end
          
           % Teleport the virus to a random location to avoid getting trapped
           virusLocation(virus, :) = PlaceCell(areaSize);
       end
   end
   % Update each immune cell in the simulation
   for immune = 1:ImmuneNum
       closest = knnsearch(virusLocation, immuneLocation(immune,:)); % Determine the closest virus to each immune cell; knnsearch is a builtin function in MATLAB
       distanceToClosest = ComputeDistance(immuneLocation(immune,:), virusLocation(closest, :));
       angleToClosest = ComputeAngle(immuneLocation(immune,:), virusLocation(closest, :)); 
       immuneLocation(immune, :) = immuneLocation(immune,:) + [cos(angleToClosest), sin(angleToClosest)] * immuneSpeed;
      
       % If immune cell captures virus, it "eats" it
       if distanceToClosest <= 2
           virusNum = virusNum - 1; % Remove virus that was eaten
           virusLocation(closest, :) = [];
           virusTarget(closest) = [];
           immuneVirusesConsumed(immune) = immuneVirusesConsumed(immune) + 1;
           % If the immune cell eats enough viruses, call in another immune cell
           if mod(immuneVirusesConsumed(immune), immuneReplicateNum) == 0
               % mod is remainder after division; if eats any multiple of 20, output will be 0 and will call in new immune cell
               immuneLocation = [immuneLocation; PlaceCell(areaSize)];
               ImmuneNum = ImmuneNum +1;
               immuneVirusesConsumed = [immuneVirusesConsumed;1];
           end
       end
   end
  
   pause(0.05) % Pause gives a delay in updating our plot so we can see each time step
   % Plot the new locations and state of our system
   plot(CellLocation(:, 1), CellLocation(:, 2), 'ok')
   hold on
   plot(ImmunoCellLocation(:, 1), ImmunoCellLocation(:, 2), 'ob')
   plot(virusLocation(:, 1), virusLocation(:,2), '.r')
   plot(immuneLocation(:, 1), immuneLocation(:,2), 'og')
   xlim([-areaSize, areaSize])
   ylim([-areaSize, areaSize])
   hold off
   % Stop looping through time if there is only one virus left
   if virusNum <= 1
       break
   end
end
%% Function to Place Agents Randomly within the Simulation Area
function location = PlaceCell(areaSize)
   location = [rand(),rand()] * 2 * areaSize - areaSize;
end
%% Function to compute the distance between two agents
function distance = ComputeDistance(point1, point2)
  % Point 1 has two values, x and y for object 1
  % Point 2 has two values, x and y for object 2
  distance = sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2);
end
%% Function to compute the angle of the line connecting the two agents (needed to update trajectory of moving agent)
function angle = ComputeAngle(point1, point2)
   angle = atan2((point2(2) - point1(2)), (point2(1) - point1(1)));
end


